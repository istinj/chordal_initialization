geometry_helpers_3d
total_least_squares_indices

function v=chordFlatten(T)
  v=zeros(12,1);
  v(1:3)=T(1,1:3)';
  v(4:6)=T(2,1:3)';
  v(7:9)=T(3,1:3)';
  v(10:12)=T(1:3,4);
end

function T=chordUnflatten(v)
  T=eye(4);
  T(1,1:3)=v(1:3)';
  T(2,1:3)=v(4:6)';
  T(3,1:3)=v(7:9)';
  T(1:3,4)=v(10:12);
end

%computes the product matrix
% such that A*B=mprod(flatten(b))*flatten(a)
% this is your jacobian vinia!
function M=mprod(v)
  M=eye(12);
  T=chordUnflatten(v);
  R=T(1:3,1:3);
  t=T(1:3,4);
  M(1:3,1:3)=R';
  M(4:6,4:6)=R';
  M(7:9,7:9)=R'; % R' on the diagonal blocks
  M(10,1:3)=t';
  M(11,4:6)=t';
  M(12,7:9)=t'; % t' on the three rows of lower part
end



% error and jacobian of a measured pose, all poses are in world frame
% input:
%   Xi: the observing robot pose (4x4 homogeneous matrix)
%   Xj: the observed robot pose (4x4 homogeneous matrix)
%   Z:   the relative transform measured between Xr1 and Xr2
%   e: 12x1 is the difference between prediction, and measurement, vectorized
%   Ji : 12x12 derivative w.r.t a the error and a perturbation of the first pose [rotation only]
%   Jj : 12x12 derivative w.r.t a the error and a perturbation of the second pose [rotation only]

function [e,Ji,Jj]=poseErrorAndJacobian(Xi,Xj,Z,i_idx,j_idx)
  Jj = eye(12,12);
  Jj = -Jj;
  Ji = mprod(chordFlatten(Z));

  % e = zeros(12,1);
  e = mprod(chordFlatten(Z))*chordFlatten(Xi)-chordFlatten(Xj);

  %ia if the pose is the first [aka the fixed one], we add a prior, saying that its
  %ia rotation matrix is the identity and it is located in the origin
  if (0==i_idx)

  end
  if (0==j_idx)
  end

 end

%linearizes the robot-robot measurements
% inputs:
%   XR: the initial robot poses (4x4xnum_poses: array of homogeneous matrices)
%   XL: the initial landmark estimates (3xnum_landmarks matrix of landmarks)
%   ZR: the robot_robot measuremenrs (4x4xnum_measurements: array of homogeneous matrices)
%   associations: 2xnum_measurements.
%                 associations(:,k)=[i_idx, j_idx]' means the kth measurement
%                 refers to an observation made from pose i_idx, that
%                 observed the pose j_idx
%   num_poses: number of poses in XR (added for consistency)
%   num_landmarks: number of landmarks in XL (added for consistency)
%   kernel_threshod: robust kernel threshold
% outputs:
%   H: the H matrix, filled
%   b: the b vector, filled
%   chi_tot: the total chi2 of the current round
%   num_inliers: number of measurements whose error is below kernel_threshold

function [H,b]=initializePoses(XR,Zr,associations,num_poses)
  global pose_dim;
  global landmark_dim;
  global flat_rotation_dimension;

  system_size=flat_rotation_dimension*num_poses;
  H=zeros(system_size, system_size);
  b=zeros(system_size,1);
  for (measurement_num=1:size(Zr,3))
    Omega=eye(12);
    Omega_i=Omega;
    Omega_j=Omega;
    Omega_ij=Omega;
    % Omega(1:9,1:9)*=1e3; % we need to pimp the rotation  part a little
    pose_i_index=associations(1,measurement_num);
    pose_j_index=associations(2,measurement_num);

    pose_i_matrix_index=chordalMatrixIndex(pose_i_index, num_poses);
    pose_j_matrix_index=chordalMatrixIndex(pose_j_index, num_poses);

    Z=Zr(:,:,measurement_num);
    Xi=XR(:,:,pose_i_index);
    Xj=XR(:,:,pose_j_index);
    [e,Ji,Jj] = poseErrorAndJacobian(Xi, Xj, Z, pose_i_matrix_index, pose_j_matrix_index);

    Jii = zeros(12,12);
    % if (1==pose_i_matrix_index)
    %   Omega_i = Omega_i*1e10;
    %   Jii = -eye(12,12);
    % end

    Jjj = zeros(12,12);
    % if (1==pose_j_matrix_index)
    %   Omega_j = Omega_j*1e10;
    %   Jjj = -eye(12,12);
    % end


    H(pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1,...
      pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1)=...
      H(pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1,...
      pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1) + Ji'*Omega_i*Ji + Jii'*Omega_i*Jii;

    H(pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1,...
      pose_j_matrix_index:pose_j_matrix_index+flat_rotation_dimension-1)=...
      H(pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1,...
      pose_j_matrix_index:pose_j_matrix_index+flat_rotation_dimension-1)+Ji'*Omega_ij*Jj;

    H(pose_j_matrix_index:pose_j_matrix_index+flat_rotation_dimension-1,...
      pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1)=...
      H(pose_j_matrix_index:pose_j_matrix_index+flat_rotation_dimension-1,...
      pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1)+Jj'*Omega_ij*Ji;

    H(pose_j_matrix_index:pose_j_matrix_index+flat_rotation_dimension-1,...
      pose_j_matrix_index:pose_j_matrix_index+flat_rotation_dimension-1)=...
      H(pose_j_matrix_index:pose_j_matrix_index+flat_rotation_dimension-1,...
      pose_j_matrix_index:pose_j_matrix_index+flat_rotation_dimension-1)+Jj'*Omega_j*Jj + Jjj'*Omega_j*Jjj;

    b(pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1)=...
      b(pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1)+Ji'*Omega_i*e;
    b(pose_j_matrix_index:pose_j_matrix_index+flat_rotation_dimension-1)=...
      b(pose_j_matrix_index:pose_j_matrix_index+flat_rotation_dimension-1)+Jj'*Omega_j*e;

  end
end

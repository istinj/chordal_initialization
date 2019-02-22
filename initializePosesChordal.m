

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

function [H,b]=initializePosesChordal(XR,Zr,associations,num_poses)
  global pose_dim;
  global landmark_dim;
  global flat_rotation_dimension;

  system_size=flat_rotation_dimension*num_poses;
  H=zeros(system_size, system_size);
  b=zeros(system_size,1);
  for (measurement_num=1:size(Zr,3))
    Omega=eye(12);
    pose_i_index=associations(1,measurement_num);
    pose_j_index=associations(2,measurement_num);

    pose_i_matrix_index=chordalMatrixIndex(pose_i_index, num_poses);
    pose_j_matrix_index=chordalMatrixIndex(pose_j_index, num_poses);

    Z=Zr(:,:,measurement_num);
    Xi=XR(:,:,pose_i_index);
    Xj=XR(:,:,pose_j_index);
    [e,Ji,Jj] = chordalErrorAndJacobian(Xi, Xj, Z);


    H(pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1,...
      pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1)=...
      H(pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1,...
      pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1) + Ji'*Omega_i*Ji;

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
      pose_j_matrix_index:pose_j_matrix_index+flat_rotation_dimension-1)+Jj'*Omega_j*Jj;

    b(pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1)=...
      b(pose_i_matrix_index:pose_i_matrix_index+flat_rotation_dimension-1)+Ji'*Omega_i*e;
    b(pose_j_matrix_index:pose_j_matrix_index+flat_rotation_dimension-1)=...
      b(pose_j_matrix_index:pose_j_matrix_index+flat_rotation_dimension-1)+Jj'*Omega_j*e;

  end    
end

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

function [H,b, chi_tot, num_inliers]=linearizePoses(XR, XL, Zr, associations,num_poses, num_landmarks, kernel_threshold)
  global pose_dim;
  global landmark_dim;
  system_size=pose_dim*num_poses+landmark_dim*num_landmarks; 
  H=zeros(system_size, system_size);
  b=zeros(system_size,1);
  chi_tot=0;
  num_inliers=0;
  for (measurement_num=1:size(Zr,3))
    Omega=eye(12);
    Omega(1:9,1:9)=Omega(1:9,1:9)*1e3; % we need to pimp the rotation  part a little
    pose_i_index=associations(1,measurement_num);
    pose_j_index=associations(2,measurement_num);
    Z=Zr(:,:,measurement_num);
    Xi=XR(:,:,pose_i_index);
    Xj=XR(:,:,pose_j_index);
    [e,Ji,Jj] = poseErrorAndJacobian(Xi, Xj, Z);
    chi=e'*Omega*e;
    if (chi>kernel_threshold)
      Omega=Omega*sqrt(kernel_threshold/chi);
      chi=kernel_threshold;
    else
      num_inliers=num_inliers+1;
    end
    chi_tot=chi_tot+chi;

    pose_i_matrix_index=poseMatrixIndex(pose_i_index, num_poses, num_landmarks);
    pose_j_matrix_index=poseMatrixIndex(pose_j_index, num_poses, num_landmarks);
    
    H(pose_i_matrix_index:pose_i_matrix_index+pose_dim-1,...
      pose_i_matrix_index:pose_i_matrix_index+pose_dim-1)=...
      H(pose_i_matrix_index:pose_i_matrix_index+pose_dim-1,...
      pose_i_matrix_index:pose_i_matrix_index+pose_dim-1)+Ji'*Omega*Ji;

    H(pose_i_matrix_index:pose_i_matrix_index+pose_dim-1,...
      pose_j_matrix_index:pose_j_matrix_index+pose_dim-1)=...
      H(pose_i_matrix_index:pose_i_matrix_index+pose_dim-1,...
      pose_j_matrix_index:pose_j_matrix_index+pose_dim-1)+Ji'*Omega*Jj;

    H(pose_j_matrix_index:pose_j_matrix_index+pose_dim-1,...
      pose_i_matrix_index:pose_i_matrix_index+pose_dim-1)=...
      H(pose_j_matrix_index:pose_j_matrix_index+pose_dim-1,...
      pose_i_matrix_index:pose_i_matrix_index+pose_dim-1)+Jj'*Omega*Ji;

    H(pose_j_matrix_index:pose_j_matrix_index+pose_dim-1,...
      pose_j_matrix_index:pose_j_matrix_index+pose_dim-1)=...
      H(pose_j_matrix_index:pose_j_matrix_index+pose_dim-1,...
      pose_j_matrix_index:pose_j_matrix_index+pose_dim-1)+Jj'*Omega*Jj;

    b(pose_i_matrix_index:pose_i_matrix_index+pose_dim-1)=...
      b(pose_i_matrix_index:pose_i_matrix_index+pose_dim-1)+Ji'*Omega*e;
    b(pose_j_matrix_index:pose_j_matrix_index+pose_dim-1)=...
      b(pose_j_matrix_index:pose_j_matrix_index+pose_dim-1)+Jj'*Omega*e;
  end
end
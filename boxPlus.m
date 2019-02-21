% implementation of the boxplus
% applies a perturbation to a set of landmarks and robot poses
% input:
%   XR: the robot poses (4x4xnum_poses: array of homogeneous matrices)
%   XL: the landmark pose (3xnum_landmarks matrix of landmarks)
%   num_poses: number of poses in XR (added for consistency)
%   num_landmarks: number of landmarks in XL (added for consistency)
%   dx: the perturbation vector of appropriate dimensions
%       the poses come first, then the landmarks
% output:
%   XR: the robot poses obtained by applying the perturbation
%   XL: the landmarks obtained by applying the perturbation

function [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx)
global pose_dim;
global landmark_dim;
for(pose_index=1:num_poses)
  pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
  dxr=dx(pose_matrix_index:pose_matrix_index+pose_dim-1);
  XR(:,:,pose_index)=v2t(dxr)*XR(:,:,pose_index);
end
for(landmark_index=1:num_landmarks)
  landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);
  dxl=dx(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,:);
  XL(:,landmark_index)=XL(:,landmark_index)+dxl;
end
end
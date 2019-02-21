

% implementation of the optimization loop with robust kernel
% applies a perturbation to a set of landmarks and robot poses
% input:
%   XR: the initial robot poses (4x4xnum_poses: array of homogeneous matrices)
%   XL: the initial landmark estimates (3xnum_landmarks matrix of landmarks)
%   Z:  the measurements (3xnum_measurements)
%   associations: 2xnum_measurements.
%                 associations(:,k)=[p_idx,l_idx]' means the kth measurement
%                 refers to an observation made from pose p_idx, that
%                 observed landmark l_idx
%   num_poses: number of poses in XR (added for consistency)
%   num_landmarks: number of landmarks in XL (added for consistency)
%   num_iterations: the number of iterations of least squares
%   damping:      damping factor (in case system not spd)
%   kernel_threshod: robust kernel threshold

% output:
%   XR: the robot poses after optimization
%   XL: the landmarks after optimization
%   chi_stats_{l,p,r}: array 1:num_iterations, containing evolution of chi2 for landmarks, projections and poses
%   num_inliers{l,p,r}: array 1:num_iterations, containing evolution of inliers landmarks, projections and poses

function [XR, XL, chi_stats_l, num_inliers_l, chi_stats_p, num_inliers_p,chi_stats_r, num_inliers_r, H, b] = doTotalLS(XR, XL, ...
  Zl, landmark_associations,...
  Zp, projection_associations,...
  Zr, pose_associations,...
  num_poses,...
  num_landmarks,...
  num_iterations,...
  damping,...
  kernel_threshold)

global pose_dim;
global landmark_dim;

chi_stats_l=zeros(1,num_iterations);
num_inliers_l=zeros(1,num_iterations);
chi_stats_p=zeros(1,num_iterations);
num_inliers_p=zeros(1,num_iterations);
chi_stats_r=zeros(1,num_iterations);
num_inliers_r=zeros(1,num_iterations);

% size of the linear system
system_size=pose_dim*num_poses+landmark_dim*num_landmarks;
for (iteration=1:num_iterations)
  H=zeros(system_size, system_size);
  b=zeros(system_size,1);
  
  if (num_landmarks)
    [H_landmarks, b_landmarks, chi_, num_inliers_] = linearizeLandmarks(XR, XL, Zl, landmark_associations,num_poses, num_landmarks, kernel_threshold);
    chi_stats_l(iteration)=chi_;
    num_inliers_l(iteration)=num_inliers_;
    
    [H_projections, b_projections, chi_, num_inliers_] = linearizeProjections(XR, XL, Zp, projection_associations,num_poses, num_landmarks, kernel_threshold);
    chi_stats_p(iteration)=chi_stats_p(iteration)+chi_;
    num_inliers_p(iteration)=num_inliers_;
  end
  
  [H_poses, b_poses, chi_, num_inliers_] = linearizePoses(XR, XL, Zr, pose_associations,num_poses, num_landmarks, kernel_threshold);
  chi_stats_r(iteration) = chi_stats_r(iteration) + chi_;
  num_inliers_r(iteration)=num_inliers_;
  
  
  H=H_poses;
  b=b_poses;
  if (num_landmarks)
    H=H+H_landmarks+H_projections;
    b=b+b_landmarks+b_projections;
  end
  
  H=H+eye(system_size)*damping;
  dx=zeros(system_size,1);
  
  % we solve the linear system, blocking the first pose
  % this corresponds to "remove" from H and b the locks
  % of the 1st pose, while solving the system
  
  dx(pose_dim+1:end)=-(H(pose_dim+1:end,pose_dim+1:end)\b(pose_dim+1:end,1));
  [XR, XL]=boxPlus(XR,XL,num_poses, num_landmarks, dx);
  
  
end
end
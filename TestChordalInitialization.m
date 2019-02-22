close all;
clear;
clc;

% This is an integrated example that comprises
%
addpath("./tools");
addpath("./tools/utilities");
addpath("./tools/visualization");
addpath("./tools/g2o_wrapper");
total_least_squares
total_least_squares_chordal_initialization

% synthesis of the virtual world
num_landmarks=100;
num_poses=10;
world_size=10;

% landmarks in a matrix, one per column
P_world=(rand(landmark_dim, num_landmarks)-0.5)*world_size;


% poses in an array of 4x4 homogeneous transform matrices
XR_true=zeros(4,4,num_poses);

% initialize 1st pose
XR_true(:,:,1)=eye(4);

% scaling coefficient for uniform random pose generation
% adjusts the translation to cover world_size
% adjusts the rotation to span the three angles;
rand_scale=eye(6);
rand_scale(1:3,1:3) = rand_scale(1:3,1:3) * (0.5*world_size);
rand_scale(4:6,4:6) = rand_scale(4:6,4:6)*pi;

for (pose_num=2:num_poses)
    xr=rand(6,1)-0.5;
    Xr=v2t(rand_scale*xr);
    XR_true(:,:,pose_num)=Xr;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POSE MEASUREMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate an odometry trajectory for the robot
num_pose_measurements=num_poses-1;
Zr=zeros(4,4,num_pose_measurements);
pose_associations=zeros(2,num_pose_measurements);

measurement_num=1;
for (pose_num=1:num_poses-1)

    Xi=XR_true(:,:,pose_num);
    Xj=XR_true(:,:,pose_num+1);
    pose_associations(:,measurement_num)=[pose_num, pose_num+1]';
    Zr(:,:,measurement_num)=inv(Xi)*Xj;

    measurement_num = measurement_num+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATION OF (WRONG) INITIAL GUESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply a perturbation to each ideal pose (construct the estimation problem)
pert_deviation=1;
pert_scale=eye(6)*pert_deviation;
XR_guess=XR_true;

for (pose_num=2:num_poses)
    xr=rand(6,1)-0.5;
    dXr=v2t(pert_scale*xr);
    XR_guess(:,:,pose_num)=dXr*XR_guess(:,:,pose_num);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALL SOLVER  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% uncomment the following to suppress pose-landmark measurements
%Zl=zeros(3,0);

% uncomment the following to suppress pose-landmark-projection measurements
%num_landmarks=0;
% Zp=zeros(3,0);

% uncomment the following to suppress pose-pose measurements
% Zr=zeros(4,4,0);
global flat_rotation_dimension;

system_size=flat_rotation_dimension*num_poses;
num_iterations = 2;
damping = 0.01;

b_vec = zeros(num_iterations, system_size);
dx_vec = zeros(num_iterations, system_size);

XR = XR_guess;

for (iteration=1:num_iterations)
  iteration
  H=zeros(system_size, system_size);
  b=zeros(system_size,1);

  [H_poses, b_poses] = initializePosesChordal(XR, Zr, pose_associations,num_poses);

  H=H_poses;
  b=b_poses;

  % H+=eye(system_size)*damping;
  dx=zeros(system_size,1);

  % we solve the linear system, blocking the first pose
  % this corresponds to "remove" from H and b the locks
  % of the 1st pose, while solving the system
%   dx(flat_rotation_dimension+1:end)=-(H(flat_rotation_dimension+1:end,flat_rotation_dimension+1:end)\b(flat_rotation_dimension+1:end,1));

  dx=H\(-b);
  XR=chordalBoxPlus(XR,num_poses, dx);
  b_vec(iteration, :) = b';
  dx_vec(iteration, :) = dx';

  % b
  % dx

end


% Plot State
figure(1);
hold on;
grid;

subplot(1,2,1);
title("Poses Initial Guess");
plot3(XR_true(1,:),XR_true(2,:),XR_true(3,:),'b*',"linewidth",2);
hold on;
plot3(XR_guess(1,:),XR_guess(2,:),XR_guess(3,:),'ro',"linewidth",2);
legend("Poses True", "Guess");grid; axis equal


subplot(1,2,2);
title("Poses After Optimization");
plot3(XR_true(1,:),XR_true(2,:),XR_true(3,:),'b*',"linewidth",2);
hold on;
plot3(XR(1,:),XR(2,:),XR(3,:),'ro',"linewidth",2);
legend("Poses True", "Guess"); grid;axis equal

%
% figure(2);
% hold on;
% grid;
% title("chi evolution");
%
% subplot(3,2,1);
% plot(chi_stats_r, 'r-', "linewidth", 2);
% legend("Chi Poses"); grid; xlabel("iterations");
% subplot(3,2,2);
% plot(num_inliers_r, 'b-', "linewidth", 2);
% legend("%inliers"); grid; xlabel("iterations");
%
% figure(3);
% title("H matrix");
% H_ =  H./H;                      % NaN and 1 element
% H_(isnan(H_))=0;                 % Nan to Zero
% H_ = abs(ones(size(H_)) - H_);   % switch zero and one
% H_ = flipud(H_);                 % switch rows
% colormap(gray(64));
% hold on;
% image([0.5, size(H_,2)-0.5], [0.5, size(H_,1)-0.5], H_*64);
% hold off;

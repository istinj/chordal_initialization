close all;
clear;
clc;

%ia add things
addpath("./tools");
addpath("./tools/utilities");
addpath("./tools/visualization");
addpath("./tools/g2o_wrapper");

global pose_dim;pose_dim=6;
global Rx0;Rx0=[0 0 0;0 0 -1;0 1 0];
global Ry0;Ry0=[0 0 1;0 0 0;-1 0 0];
global Rz0;Rz0=[0 -1 0;1 0 0;0 0 0];
global flat_rotation_dimension;flat_rotation_dimension = 9;

color_cyan = [0 .75 .75];
color_red = [.9 .1 .1];
color_green = [.1 .9 .1];
color_grey = [.1 .1 .1];

%ia ground truth generation
num_poses = 200;
world_size = 50;

num_iterations = 1

zero_guess = true
meas_have_noise = true

meas_noise_sigma_trans = 1.0;
meas_noise_sigma_rot = 0.1;

meas_noise_sigma = eye(pose_dim,pose_dim);
meas_noise_sigma(1:3,1:3) = meas_noise_sigma(1:3,1:3)*meas_noise_sigma_trans;
meas_noise_sigma(4:end,4:end) = meas_noise_sigma(4:end,4:end)*meas_noise_sigma_rot;
meas_noise_sigma

Factor = struct("id",-1,"pose",eye(flat_rotation_dimension,flat_rotation_dimension),"ass",[-1;-1],...
  "hii_idx_start",-1,"hii_idx_end",-1,...
  "hjj_idx_start",-1,"hjj_idx_end",-1);

%% Data generation

Xr_gt = zeros(4,4,num_poses);
Xr_gt(:,:,1) = eye(4,4); %ia first pose is in the origin

%ia random poses in the world
rand_scale=eye(6);
rand_scale(1:3,1:3) = rand_scale(1:3,1:3) * (0.5*world_size);
rand_scale(4:6,4:6) = rand_scale(4:6,4:6)*pi;
for p=2:num_poses
  xr=rand(6,1)-0.5;
  Xr=v2t(rand_scale*xr);
  Xr_gt(:,:,p)=Xr;
end

%ia generate pose measurements and data association
num_meas = num_poses-1;
factors = {num_poses};
for m=1:num_meas
  Xr_from = Xr_gt(:,:,m);
  Xr_to = Xr_gt(:,:,m+1);
  
  factors{m} = Factor();
  factors{m}.id = m;
  factors{m}.pose = inv(Xr_from)*Xr_to;
  factors{m}.ass = [m;m+1];
  factors{m}.hii_idx_start = chordalMatrixIndex(m, num_poses);
  factors{m}.hii_idx_end = chordalMatrixIndex(m, num_poses)+flat_rotation_dimension-1;
  factors{m}.hjj_idx_start = chordalMatrixIndex(m+1, num_poses);
  factors{m}.hjj_idx_end = chordalMatrixIndex(m+1, num_poses)+flat_rotation_dimension-1;
  
  if meas_have_noise
    R = chol(meas_noise_sigma);
    dz_noise = (zeros(1,pose_dim) + randn(1,pose_dim)*R)';
    disp('noise = ');
    dz_noise'
    factors{m}.pose = v2t(dz_noise)*factors{m}.pose;
  end
end

%ia generate a wrong initial guess
Xr_guess = zeros(4,4,num_poses);
Xr_guess(:,:,1) = Xr_gt(:,:,1); %ia fix the first pose
pert_deviation=world_size;
pert_scale=eye(6)*pert_deviation;
for p=2:num_poses
  if zero_guess
    %ia zero guess;
    Xr_guess(:,:,p) = eye(4,4);
  else
    %ia random guess
    xr=rand(6,1)-0.5;
    dxr = pert_scale*xr;
    Xr_guess(:,:,p) = v2t(dxr)*Xr_gt(:,:,p);
  end
end


%% Solve Rotations [positions are fixed]
%ia start doing things
Xr = Xr_guess;

for iteration=1:num_iterations
  H = zeros(flat_rotation_dimension*num_poses,flat_rotation_dimension*num_poses);
  b = zeros(flat_rotation_dimension*num_poses,1);
  
  Omega = eye(flat_rotation_dimension,flat_rotation_dimension);
  % Omega(1:9,1:9) = Omega(1:9,1:9)*1e9;
  % Omega(10:end,10:end) = zeros(3,3);
  % Omega;
  for z=1:num_meas
    factor = factors{z};
    Ji = zeros(flat_rotation_dimension,flat_rotation_dimension);
    Jj = zeros(flat_rotation_dimension,flat_rotation_dimension);
    Z = factor.pose;
    Xi = Xr(:,:,factor.ass(1));
    Xj = Xr(:,:,factor.ass(2));
    
    flat_xi = flattenByRow(Xi);
    flat_xj = flattenByRow(Xj);
    
    flat_xi = flat_xi(1:9);
    flat_xj = flat_xj(1:9);
    
    M = zeros(9,9);
    M(1:3,1:3) = Z(1:3,1:3)';
    M(4:6,4:6) = Z(1:3,1:3)';
    M(7:9,7:9) = Z(1:3,1:3)';
    
    %ia error: Zij - hij(X) = Zij - inv(Xi)Xj = XiZij - Xj -> applying magic
    %ia       -> Mij*flat(Xi) - flat(Xj)
    %ia jacobians: Ji = Mij; Jj = -I(flat_rotation_dimension x flat_rotation_dimension);
    error = zeros(9,1);
    error = M * flat_xi - flat_xj;
    Ji = M;
    Jj = -1.0 * eye(flat_rotation_dimension,flat_rotation_dimension);
    
    %ia contruct least squares
    %ia diagonal components
    H(factor.hii_idx_start:factor.hii_idx_end,factor.hii_idx_start:factor.hii_idx_end) = ...
      H(factor.hii_idx_start:factor.hii_idx_end,factor.hii_idx_start:factor.hii_idx_end)+Ji'*Omega*Ji;
    H(factor.hjj_idx_start:factor.hjj_idx_end,factor.hjj_idx_start:factor.hjj_idx_end) = ...
      H(factor.hjj_idx_start:factor.hjj_idx_end,factor.hjj_idx_start:factor.hjj_idx_end)+Jj'*Omega*Jj;
    %ia mixed components
    H(factor.hii_idx_start:factor.hii_idx_end,factor.hjj_idx_start:factor.hjj_idx_end) = ...
      H(factor.hii_idx_start:factor.hii_idx_end,factor.hjj_idx_start:factor.hjj_idx_end)+Ji'*Omega*Jj;
    H(factor.hjj_idx_start:factor.hjj_idx_end,factor.hii_idx_start:factor.hii_idx_end) = ...
      H(factor.hjj_idx_start:factor.hjj_idx_end,factor.hii_idx_start:factor.hii_idx_end)+Jj'*Omega*Ji;
    %ia rhs vector
    b(factor.hii_idx_start:factor.hii_idx_end,1) = ...
      b(factor.hii_idx_start:factor.hii_idx_end,1)+Ji'*Omega*error;
    b(factor.hjj_idx_start:factor.hjj_idx_end,1) = ...
      b(factor.hjj_idx_start:factor.hjj_idx_end,1)+Jj'*Omega*error;
    
  end
  
  %   H(1:flat_rotation_dimension,1:flat_rotation_dimension) = ...
  %     H(1:flat_rotation_dimension,1:flat_rotation_dimension)+eye(flat_rotation_dimension, flat_rotation_dimension);
  %
  %   b(1:flat_rotation_dimension) = b(1:flat_rotation_dimension)+flattenByRow(Xr(:,:,1)-Xr_gt(:,:,1));
  
  %ia compute solution of the system
  dx = zeros(flat_rotation_dimension*num_poses,1);
  %   dx = H\(-b); %ia simple solution
  dx(flat_rotation_dimension+1:end)=-(H(flat_rotation_dimension+1:end,flat_rotation_dimension+1:end)\b(flat_rotation_dimension+1:end,1));%ia discard first pose
  
  b;
  dx;
  rank(H)
  
  %ia apply increment
  for p=1:num_poses
    idx_start = chordalMatrixIndex(p,num_poses);
    idx_end = idx_start+flat_rotation_dimension-1;
    increment = dx(idx_start:idx_end,1);
    
    xr = flattenByRow(Xr(:,:,p));
    xr(1:flat_rotation_dimension,1) = xr(1:flat_rotation_dimension,1)+increment;
    Xr(:,:,p) = unflattenByRow(xr, true); %ia second argument re-enforces orthogonality constraint
  end
  
  %   pause;
end


%% Solve Full Pose [with Chordal?]
for iteration=1:num_iterations
  H = zeros(pose_dim*num_poses,pose_dim*num_poses);
  b = zeros(pose_dim*num_poses,1);
  dx = zeros(pose_dim*num_poses,1);
  Omega = eye(12,12);
  
  chi_tot = 0;
  kernel_threshold = 1e10;
  num_inliers = 0;
  for z=1:num_meas
    factor = factors{z};
    Z = factor.pose;
    Xi = Xr(:,:,factor.ass(1));
    Xj = Xr(:,:,factor.ass(2));
    
    [e,Ji,Jj]=poseErrorAndJacobian(Xi,Xj,Z);
    chi=e'*Omega*e;
    if (chi>kernel_threshold)
      Omega=Omega*sqrt(kernel_threshold/chi);
      chi=kernel_threshold;
    else
      num_inliers=num_inliers+1;
    end
    chi_tot=chi_tot+chi;
    
    pose_i_matrix_index=poseMatrixIndex(factor.ass(1), num_poses, 0);
    pose_j_matrix_index=poseMatrixIndex(factor.ass(2), num_poses, 0);
    
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
  
  rank(H)
  
  %   dx = H\(-b); %ia simple
  dx(pose_dim+1:end,1) = H(pose_dim+1:end,pose_dim+1:end)\(-b(pose_dim+1:end,1)); %ia discard first pose
  
  for(p=1:num_poses)
    pose_matrix_index=poseMatrixIndex(p, num_poses, 0);
    dxr=dx(pose_matrix_index:pose_matrix_index+pose_dim-1);
    Xr(:,:,p)=v2t(dxr)*Xr(:,:,p);
  end
  
end

%% Show things
%ia show the created world
figure("Name","Chordal Initialization Dioporco");
subplot(1,2,1);
title("Ground Truth and Initial Guess");
hold on; grid; axis([-world_size*0.5, world_size*0.5,-world_size*0.5, world_size*0.5]);
px = Xr_gt(1,4,:);
py = Xr_gt(2,4,:);
pz = Xr_gt(3,4,:);
scatter3(px,py,pz,30,'filled','MarkerFaceColor',color_cyan);
px = Xr_guess(1,4,:);
py = Xr_guess(2,4,:);
pz = Xr_guess(3,4,:);
scatter3(px,py,pz,100,'o','MarkerEdgeColor',color_red);
legend("GT", "Guess");
hold off;

%ia show the updated world
subplot(1,2,2);
title("Ground Truth and Initialization");
hold on; grid; axis([-world_size*0.5, world_size*0.5,-world_size*0.5, world_size*0.5]);
px = Xr_gt(1,4,:);
py = Xr_gt(2,4,:);
pz = Xr_gt(3,4,:);
scatter3(px,py,pz,30,'filled','MarkerFaceColor',color_cyan);
px = Xr(1,4,:);
py = Xr(2,4,:);
pz = Xr(3,4,:);
scatter3(px,py,pz,100,'o','MarkerEdgeColor',color_grey);
legend("GT", "Initialization");
hold off;









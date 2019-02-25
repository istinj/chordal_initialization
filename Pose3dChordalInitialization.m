close all;
clear;
clc;

global flat_rotation_dimension;flat_rotation_dimension = 12;

chordal_dimension = 12;
zero_guess = 0

Factor = struct("id",-1,"pose",eye(chordal_dimension,chordal_dimension),"ass",[-1;-1],...
  "hii_idx_start",-1,"hii_idx_end",-1,...
  "hjj_idx_start",-1,"hjj_idx_end",-1);

%% Compute things
%ia add things
addpath("./tools");
addpath("./tools/utilities");
addpath("./tools/visualization");
addpath("./tools/g2o_wrapper");

color_cyan = [0 .75 .75];
color_red = [.9 .1 .1];
color_green = [.1 .9 .1];

%ia ground truth generation
num_poses = 200;
world_size = 20;
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
  factors{m}.hii_idx_end = chordalMatrixIndex(m, num_poses)+chordal_dimension-1;
  factors{m}.hjj_idx_start = chordalMatrixIndex(m+1, num_poses);
  factors{m}.hjj_idx_end = chordalMatrixIndex(m+1, num_poses)+chordal_dimension-1;
end

%ia generate a wrong initial guess
Xr_guess = zeros(4,4,num_poses);
Xr_guess(:,:,1) = Xr_gt(:,:,1); %ia fix the first pose
pert_deviation=10;
pert_scale=eye(6)*pert_deviation;
for p=2:num_poses
  if zero_guess
    %ia zero guess;
    Xr_guess(:,:,p) = eye(4,4);
  else
    %ia random guess
    xr=rand(6,1)-0.5;
    Xr_guess(:,:,p) = v2t(pert_scale*xr)*Xr_gt(:,:,p);
  end
end

%ia start doing things
Xr = Xr_guess;

H = zeros(chordal_dimension*num_poses,chordal_dimension*num_poses);
b = zeros(chordal_dimension*num_poses,1);

Omega = eye(chordal_dimension,chordal_dimension);
% Omega(1:9,1:9) = Omega(1:9,1:9)*1e9;
% Omega(10:end,10:end) = zeros(3,3);
% Omega;
for z=1:num_meas
  factor = factors{z};
  Ji = zeros(chordal_dimension,chordal_dimension);
  Jj = zeros(chordal_dimension,chordal_dimension);
  Z = factor.pose;
  Xi = Xr(:,:,factor.ass(1));
  Xj = Xr(:,:,factor.ass(2));
  
  flat_xi = flattenByRow(Xi);
  flat_xj = flattenByRow(Xj);
  Mij = mprod(flattenByRow(Z)); %ia creates a chordal_dimensionxchordal_dimension matrix from the isometry
  
  %ia error: Zij - hij(X) = Zij - inv(Xi)Xj = XiZij - Xj -> applying magic
  %ia       -> Mij*flat(Xi) - flat(Xj)
  %ia jacobians: Ji = Mij; Jj = -I(chordal_dimensionxchordal_dimension);
  error = Mij * flat_xi - flat_xj;
  Ji = Mij;
  Jj = -1.0 * eye(chordal_dimension,chordal_dimension);
  
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

% H(1:chordal_dimension,1:chordal_dimension) = ...
%   H(1:chordal_dimension,1:chordal_dimension)+eye(chordal_dimension, chordal_dimension);
%
% b(1:chordal_dimension) = b(1:chordal_dimension)+flattenIsometry(Xr(:,:,1)-Xr_gt(:,:,1));

%ia compute solution of the system
dx = zeros(chordal_dimension*num_poses,1);
% dx = H\(-b); %ia simple solution
dx(chordal_dimension+1:end)=-(H(chordal_dimension+1:end,chordal_dimension+1:end)\b(chordal_dimension+1:end,1));%ia discard first pose

% b
% dx
% rank(H)

%ia apply increment
for p=1:num_poses
  idx_start = chordalMatrixIndex(p,num_poses);
  idx_end = idx_start+chordal_dimension-1;
  increment = dx(idx_start:idx_end,1);
  
  %ia apply increment to translation
  Xr(1:3,4,p) = Xr(1:3,4,p) + increment(10:end);
  
  %ia apply increment to rotation
  Xr(1:3,1,p) = Xr(1:3,1,p) + increment(1:3);
  Xr(1:3,2,p) = Xr(1:3,2,p) + increment(4:6);
  Xr(1:3,3,p) = Xr(1:3,3,p) + increment(7:9);
  %ia re-enforce the rotation constraint
  [U,S,V] = svd(Xr(1:3,1:3));
  R_enforced = U*V';
  Xr(1:3,1:3) = R_enforced;
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
scatter3(px,py,pz,100,'o','MarkerEdgeColor',color_red);
legend("GT", "Initialization");
hold off;









% error and jacobian of a measured pose, all poses are in world frame
% input:
%   Xi: the observing robot pose (4x4 homogeneous matrix)
%   Xj: the observed robot pose (4x4 homogeneous matrix)
%   Z:   the relative transform measured between Xr1 and Xr2
%   e: 12x1 is the difference between prediction, and measurement, vectorized
%   Ji : 12x6 derivative w.r.t a the error and a perturbation of the
%       first pose
%   Jj : 12x6 derivative w.r.t a the error and a perturbation of the
%       second pose

function [e,Ji,Jj]=poseErrorAndJacobianGeodesic(Xi,Xj,Z)
global pose_dim;
Ji=zeros(6,6);
Jj=zeros(6,6);
e=zeros(6,1);

Zhat = inv(Xi)*Xj;
E = inv(Z)*Zhat;
e = t2v_euler(E);

%numerical jacobians?
eps = 1e-7;
for u=1:pose_dim
  delta = zeros(pose_dim,1);
  delta(pose_dim,1) = eps;
  %ia things
end

end
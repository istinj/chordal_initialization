
% error and jacobian of a measured landmark
% input:
%   Xr: the robot pose in world frame (4x4 homogeneous matrix)
%   Xl: the landmark pose (3x1 vector, 3d pose in world frame)
%   z:  projection of the landmark on the image plane
% output:
%   e: 2x1 is the difference between prediction and measurement
%   Jr: 2x6 derivative w.r.t a the error and a perturbation on the
%       pose
%   Jl: 2x3 derivative w.r.t a the error and a perturbation on the
%       landmark
%   is_valid: true if projection ok

function [is_valid, e,Jr,Jl]=projectionErrorAndJacobian(Xr,Xl,z)
  global K;
  global image_rows;
  global image_cols;
  is_valid=false;
  e=[0;0];
  Jr=zeros(2,6);
  Jl=zeros(2,3);
  
  % inverse transform
  iR=Xr(1:3,1:3)';
  it=-iR*Xr(1:3,4);

  pw=iR*Xl+it; %point prediction, in world scale
  if (pw(2)<0)
     return;
  end

  Jwr=zeros(3,6);
  Jwr(1:3,1:3)=-iR;
  Jwr(1:3,4:6)=iR*skew(Xl);
  Jwl=iR;
  
  p_cam=K*pw;
  iz=1./p_cam(3);
  z_hat=p_cam(1:2)*iz;
  if (z_hat(1)<0 || ...
      z_hat(1)>image_cols ||...
      z_hat(2)<0 || ...
      z_hat(2)>image_rows)
    return;
  end

  iz2=iz*iz;
  Jp=[iz, 0, -p_cam(1)*iz2;
      0, iz, -p_cam(2)*iz2];
  
  e=z_hat-z;
  Jr=Jp*K*Jwr;
  Jl=Jp*K*Jwl;
  is_valid=true;
end
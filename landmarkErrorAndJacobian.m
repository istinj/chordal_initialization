
% error and jacobian of a measured landmark
% input:
%   Xr: the robot pose in world frame (4x4 homogeneous matrix)
%   Xl: the landmark pose (3x1 vector, 3d pose in world frame)
%   z:  measured position of landmark
% output:
%   e: 3x1 is the difference between prediction and measurement
%   Jr: 3x6 derivative w.r.t a the error and a perturbation on the
%       pose
%   Jl: 3x3 derivative w.r.t a the error and a perturbation on the
%       landmark

function [e,Jr,Jl]=landmarkErrorAndJacobian(Xr,Xl,z)
  % inverse transform
  iR=Xr(1:3,1:3)';
  it=-iR*Xr(1:3,4);
  
  %prediction
  z_hat=iR*Xl+it; 
  e=z_hat-z;
  Jr=zeros(3,6);
  Jr(1:3,1:3)=-iR;
  Jr(1:3,4:6)=iR*skew(Xl);
  Jl=iR;
end

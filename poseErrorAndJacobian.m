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

function [e,Ji,Jj]=poseErrorAndJacobian(Xi,Xj,Z)
  global Rx0;
  global Ry0;
  global Rz0;
  Ri=Xi(1:3,1:3);
  Rj=Xj(1:3,1:3);
  ti=Xi(1:3,4);
  tj=Xj(1:3,4);
  tij=tj-ti;
  Ri_transpose=Ri';
  Ji=zeros(12,6);
  Jj=zeros(12,6);
  
  dR_dax=Ri_transpose*Rx0*Rj;
  dR_day=Ri_transpose*Ry0*Rj;
  dR_daz=Ri_transpose*Rz0*Rj;
  
  Jj(1:9,4)=reshape(dR_dax, 9, 1);
  Jj(1:9,5)=reshape(dR_day, 9, 1);
  Jj(1:9,6)=reshape(dR_daz, 9, 1);
  Jj(10:12,1:3)=Ri_transpose;
  
  Jj(10:12,4:6)=-Ri_transpose*skew(tj);
  Ji=-Jj;

  Z_hat=eye(4);
  Z_hat(1:3,1:3)=Ri_transpose*Rj;
  Z_hat(1:3,4)=Ri_transpose*tij;
  e=flattenIsometryByColumns(Z_hat-Z);
 end
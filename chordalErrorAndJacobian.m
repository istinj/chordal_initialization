
% error and jacobian of a measured pose, all poses are in world frame
% input:
%   Xi: the observing robot pose (4x4 homogeneous matrix)
%   Xj: the observed robot pose (4x4 homogeneous matrix)
%   Z:   the relative transform measured between Xr1 and Xr2
%   e: 12x1 is the difference between prediction, and measurement, vectorized
%   Ji : 12x12 derivative w.r.t a the error and a perturbation of the first pose [rotation only]
%   Jj : 12x12 derivative w.r.t a the error and a perturbation of the second pose [rotation only]

function [e,Ji,Jj]=chordalErrorAndJacobian(Xi,Xj,Z,i_idx,j_idx)
  Jj = eye(12,12);
  Jj = -Jj;
  Ji = mprod(flattenIsometry(Z));

  % e = zeros(12,1);
  e = mprod(flattenIsometry(Z))*flattenIsometry(Xi)-flattenIsometry(Xj);

  %ia if the pose is the first [aka the fixed one], we add a prior, saying that its
  %ia rotation matrix is the identity and it is located in the origin
  if (0==i_idx)

  end
  if (0==j_idx)
  end

 end
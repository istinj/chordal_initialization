
function XR=chordalBoxPlus(XR, num_poses, dx)
global pose_dim;
global landmark_dim;
global flat_rotation_dimension;
for(pose_index=1:num_poses)
  pose_matrix_index=chordalMatrixIndex(pose_index, num_poses);
  dxr=dx(pose_matrix_index:pose_matrix_index+flat_rotation_dimension-1);
  XR(1:3,1) = XR(1:3,1)+dxr(1:3);
  XR(1:3,2) = XR(1:3,2)+dxr(4:6);
  XR(1:3,3) = XR(1:3,3)+dxr(7:9);
  
  [U,S,V] = svd(XR(1:3,1:3));
  R_enforced = U*V';
  XR(1:3,1:3) = R_enforced;
  
  XR(1:3,4) = XR(1:3,4) + dxr(10:end);
end
end
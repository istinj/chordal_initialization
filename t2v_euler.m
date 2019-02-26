function v=t2v(T)
global pose_dim;
eulXYZ = rotm2eul(T(1:3,1:3),'XYZ');
v = zeros(pose_dim,1);
v(1:3,1) = T(1:3,4);
v(4:end,1) = eulXYZ';
end
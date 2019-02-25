function X_spiral = createSpiral(radius,horizontal_steps,vertical_steps)

num_poses = horizontal_steps*vertical_steps;
X_spiral = zeros(4,4,num_poses);
X_spiral(:,:,1) = eye(4,4);

p = 1;
for h = 1:horizontal_steps
  for v = 1:vertical_steps
    X = eye(4,4);
    
    rotz = (-pi + 2*v*pi / vertical_steps) * [0 0 1];
    roty = -0.5*pi + p*pi / (horizontal_steps*vertical_steps) * [0 1 0];
    p = p+1;
    
    Rz = rotationVectorToMatrix(rotz);
    Ry = rotationVectorToMatrix(roty);
    R = Rz*Ry;
    t = R*[radius;0;0];
    X(1:3,1:3) = R;
    X(1:3,4) = t;
    
    X_spiral(:,:,p) = X
  end
end
end
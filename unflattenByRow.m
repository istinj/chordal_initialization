function T=unflattenByRow(v,enforce_rotations)
  T=eye(4);
  T(1,1:3)=v(1:3)';
  T(2,1:3)=v(4:6)';
  T(3,1:3)=v(7:9)';
  T(1:3,4)=v(10:12);
  
  if enforce_rotations
    [U,S,V] = svd(T(1:3,1:3));
    R_enforced = U*V';
    T(1:3,1:3) = R_enforced;
  end
end
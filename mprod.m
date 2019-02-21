%computes the product matrix
% such that A*B=mprod(flatten(b))*flatten(a)
% this is your jacobian vinia!
function M=mprod(v)
  M=eye(12);
  T=unflattenIsometry(v);
  R=T(1:3,1:3);
  t=T(1:3,4);
  M(1:3,1:3)=R';
  M(4:6,4:6)=R';
  M(7:9,7:9)=R'; % R' on the diagonal blocks
  M(10,1:3)=t';
  M(11,4:6)=t';
  M(12,7:9)=t'; % t' on the three rows of lower part
end
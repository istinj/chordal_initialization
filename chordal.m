function v=flatten(T)
  v=zeros(12,1);
  v(1:3)=T(1,1:3)';
  v(4:6)=T(2,1:3)';
  v(7:9)=T(3,1:3)';
  v(10:12)=T(1:3,4);
end;

function T=unflatten(v)
  T=eye(4);
  T(1,1:3)=v(1:3)';
  T(2,1:3)=v(4:6)';
  T(3,1:3)=v(7:9)';
  T(1:3,4)=v(10:12);
end;

%computes the product matrix
% such that A*B=mprod(flatten(b))*flatten(a)
% this is your jacobian vinia!
function M=mprod(v)
  M=eye(12);
  T=unflatten(v)
  R=T(1:3,1:3)
  t=T(1:3,4);
  M(1:3,1:3)=M(4:6,4:6)=M(7:9,7:9)=R'; % R' on the diagonal blocks
  M(10,1:3)=M(11,4:6)=M(12,7:9)=t';    % t' on the three rows of lower part
end

function doTest()
  v1=rand(12,1)
  v2=rand(12,1)
  T1=unflatten(v1)
  T2=unflatten(v2)
  T12=T1*T2
  T12x=unflatten(mprod(v2)*v1)
  E=T12-T12x

  T2_inv = [T2(1:3,1:3)', -T2(1:3,1:3)'*T2(1:3,4);zeros(1,3),1]
  T1_inv = [T1(1:3,1:3)', -T1(1:3,1:3)'*T1(1:3,4);zeros(1,3),1]
  Tdiff = T1_inv*T2
  v1_inv = flatten(T1_inv)
  Tdiff_x = unflatten(mprod(v2)*v1_inv)
end

function v=flattenByCol(T)
  v=zeros(12,1);
  v(1:3)=T(1:3,1);
  v(4:6)=T(1:3,2);
  v(7:9)=T(1:3,3);
  v(10:12)=T(1:3,4);
end
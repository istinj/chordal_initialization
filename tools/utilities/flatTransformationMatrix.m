    function M=flatTransformationMatrix(v)
      T=unflattenIsometry(v);
      R=T(1:3,1:3);
      t=T(1:3,4);
      M=eye(12);
      M(1:3,1:3)=R';
      M(4:6,4:6)=R';
      M(7:9,7:9)=R';
      M(10,1:3)=t';
      M(11,4:6)=t';
      M(12,7:9)=t';
    end
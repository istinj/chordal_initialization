    function T=unflattenIsometryByColumns(v)
      T=eye(4);
      T(1:3,1:3)=reshape(v(1:9),3,3);
      T(1:3,4)=v(10:12);
    end
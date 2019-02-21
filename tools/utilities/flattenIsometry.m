    function v=flattenIsometry(T)
      v=zeros(12,1);
      v(1:9)=reshape(T(1:3,1:3)',9,1);
      v(10:12)=T(1:3,4);
    end
    %from 6d vector to homogeneous matrix
    function T=v2t(v)
      T=eye(4);
      T(1:3,1:3)=Rx(v(4))*Ry(v(5))*Rz(v(6));
      T(1:3,4)=v(1:3);
    end
    %derivative of rotation matrix around z
    function R=Rx_prime(rot_x)
      dc=-sin(rot_x); %derivative of cos(rot(x)
      ds=cos(rot_x);  %derivative of sin(rot(x)
      R= [0  0  0;
        0  dc  -ds;
        0  ds  dc];
    end
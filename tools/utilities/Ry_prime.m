    %derivative of rotation matrix around y
    function R=Ry_prime(rot_y)
      dc=-sin(rot_y); %derivative of cos(rot(x)
      ds=cos(rot_y);  %derivative of sin(rot(x)
      R= [dc  0 ds;
        0  0  0;
        -ds  0 dc];
    end
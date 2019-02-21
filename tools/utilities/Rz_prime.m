    %derivative of rotation matrix around z
    function R=Rz_prime(rot_z)
      dc=-sin(rot_z); %derivative of cos(rot(x)
      ds=cos(rot_z);  %derivative of sin(rot(x)
      R= [ dc  -ds  0;
        ds  dc  0;
        0  0  0];
    end
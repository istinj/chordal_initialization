    %rotation matrix around y axis
    function R=Ry(rot_y)
      c=cos(rot_y);
      s=sin(rot_y);
      R= [c  0  s;
        0  1  0;
        -s  0 c];
    end
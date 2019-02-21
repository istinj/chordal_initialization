%rotation matrix around x axis
function R=Rx(rot_x)
c=cos(rot_x);
s=sin(rot_x);
R= [1  0  0;
  0  c  -s;
  0  s  c];
end
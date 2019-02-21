%rotation matrix around z axis
function R=Rz(rot_z)
c=cos(rot_z);
s=sin(rot_z);
R= [ c  -s  0;
  s  c  0;
  0  0  1];
end
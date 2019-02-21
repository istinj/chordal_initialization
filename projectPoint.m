
% projects a point
function p_img=projectPoint(Xr,Xl)
  global image_cols;
  global image_rows;
  global K;
  iXr=inv(Xr);
  p_img=[-1;-1];
  pw=iXr(1:3,1:3)*Xl+iXr(1:3,4);
  if (pw(3)<0)
     return;
  end
  p_cam=K*pw;
  iz=1./p_cam(3);
  p_cam=p_cam*iz;
  if (p_cam(1)<0 || ...
      p_cam(1)>image_cols ||...
      p_cam(2)<0 || ...
      p_cam(2)>image_rows)
    return;
  end
  p_img=p_cam(1:2);
end
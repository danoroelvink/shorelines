function [dir0] = dirintd(x,hs,dir,x0)
%Interpolaton of directions
cs=interp1(x,hs.*cosd(dir),x0);
sn=interp1(x,hs.*sind(dir),x0);
dir0=atan2d(sn,cs);
end


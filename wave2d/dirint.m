function [dir0] = dirintd(x,hs,dir,x0)
%Interpolaton of directions
cs=interp1(x,hs.*cos(dir),x0);
sn=interp1(x,hs.*sin(dir),x0);
dir0=atan2(sn,cs);
end


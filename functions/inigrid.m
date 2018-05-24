function [ xg,yg,zg ] = inigrid( dx,dy,Lx,Ly, zdeep,slope,zshallow )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
x0=0:dx:Lx;
y0=0:dy:Ly;
for j=1:length(y0);
    xg(j,:)=x0;
end
for i=1:length(x0);
    yg(:,i)=y0;
end
zg=min(zdeep+xg*slope,zshallow);

end


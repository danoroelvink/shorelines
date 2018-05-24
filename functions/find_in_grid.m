function [row,col]=find_in_grid(x,y,x0,y0);
range=50000;
inrange=abs(x-x0)<range&abs(y-y0)<range;
dist=zeros(size(x));
dist(inrange)=hypot(x(inrange)-x0,y(inrange)-y0);
dist(~inrange)=1.e10;
mindist=min(min(dist));
[row,col,val]=find(dist==mindist);


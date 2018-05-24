function [xi,yi]=select_polygon
% Select polygon to include in bathy
hold on;
xi = [];yi=[];
n = 0;
% Loop, picking up the points.
but = 1;
while but == 1
    n = n+1;
    [xi(n),yi(n),but] = ginput(1);
    plot(xi,yi,'r-o');
end

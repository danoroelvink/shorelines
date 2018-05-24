function [x_mc,y_mc]=insert_section(x,y,x_mc,y_mc,i_mc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[ xold,yold,n_mc,i1,i2 ] = get_one_polygon( x_mc,y_mc,i_mc );%% insert x,y back into x_mc,y_mc
x_mc=[x_mc(1:i1-1),x,x_mc(i2+1:end)];
y_mc=[y_mc(1:i1-1),y,y_mc(i2+1:end)];

end


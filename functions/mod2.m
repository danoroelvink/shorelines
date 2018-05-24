function [ z ] = mod2( x,y )
%% modified mod function that returns y instead of 0 if x==y
z=mod(x,y);
z(z==0)=y;
end


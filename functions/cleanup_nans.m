function [x_mc,y_mc]=cleanup_nans(x_mc,y_mc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
i=1;
while i<=length(x_mc)
   if i==1&& isnan(x_mc(i))
       x_mc=x_mc(2:end);
       y_mc=y_mc(2:end);
   elseif i==length(x_mc)&& isnan(x_mc(i))
       x_mc=x_mc(1:end-1);
       y_mc=y_mc(1:end-1);
   elseif i>1&&i<length(x_mc)&&isnan(x_mc(i-1))&&isnan(x_mc(i+1))
       x_mc=[x_mc(1:i-1),x_mc(i+2:end)];
       y_mc=[y_mc(1:i-1),y_mc(i+2:end)];
   elseif i>1&&i<length(x_mc)-1&&isnan(x_mc(i-1))&&isnan(x_mc(i+2))
       x_mc=[x_mc(1:i-1),x_mc(i+3:end)];
       y_mc=[y_mc(1:i-1),y_mc(i+3:end)];
   elseif i==length(x_mc)-1&&isnan(x_mc(i))
       x_mc=[x_mc(1:i-1)];
       y_mc=[y_mc(1:i-1)];      
   elseif isnan(x_mc(i))&&isnan(x_mc(i+1))
       x_mc=[x_mc(1:i),x_mc(i+2:end)];
       y_mc=[y_mc(1:i),y_mc(i+2:end)];
   elseif i>1&&x_mc(i)==x_mc(i-1)&&y_mc(i)==y_mc(i-1)
       x_mc=[x_mc(1:i-1),x_mc(i+1:end)];
       y_mc=[y_mc(1:i-1),y_mc(i+1:end)];
       
   else
       i=i+1;
   end
end


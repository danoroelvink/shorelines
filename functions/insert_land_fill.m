function [x,y,xmax,ymax,xmin,ymin,nmax] = insert_land_fill(x,y,ld,it,xmax,ymax,xmin,ymin,i_mc,n_mc,nmax)
%function [x,y,nmax] = insert_land_fill(x,y,ld,it,i_mc,n_mc,nmax)
 %insert land area behind shoreline , the area width ld (m)
 %surrounding with black line
 % work with fill_sections function

  
 if it==0  % this value should store
     xmin(i_mc)=min(x);
     ymin(i_mc)=min(y);
     xmax(i_mc)=max(x);
     ymax(i_mc)=max(y);
     nmax=n_mc;
 end
%  if n_mc > nmax
%      nads=n_mc-nmax
%      if i_mc > nmax
%          for ina=1:nads
%              xmin(i_mc)=min(x);
%              ymin(i_mc)=min(y);
%              xmax(i_mc)=max(x);
%              ymax(i_mc)=max(y);
%          end
%      end 
%  end
       
 
 
    o=atan2((y(end)-y(1)),(x(end)-x(1)))*180/pi;
if o >= -45 &&  o<=45   %% Shoreline horizontal & land on the right side (down)
     x=[xmin(i_mc),x,xmax(i_mc)];
     y=[min(ymin)-ld,y,min(ymin)-ld];
elseif abs (o) >= 135  %% Shoreline horizontal & land on the left side (up)
     x=[xmax(i_mc),x,xmin(i_mc)];
     y=[max(ymax)+ld,y,max(ymax)+ld];
elseif o > 45 &&  o<135    %% vertical Shoreline  & land on the right side 
    x=[max(xmax)+ld,x,max(xmax)+ld];
     y=[ymin(i_mc),y,ymax(i_mc)];
elseif o < -45 &&  o>-135   %% vertical Shoreline  & land on the left side 
     x=[min(xmin)-ld,x,min(xmin)-ld];
     y=[ymax(i_mc),y,ymin(i_mc)];
end

%     horizontal = abs(xmax(i_mc)-xmin(i_mc))>abs(ymax(i_mc)-ymin(i_mc))   % the shoreline oreintation
% if horizontal && x(end)>x(1)   %% Shoreline horizontal & land on the right side (down)
%      x=[xmin(i_mc),x,xmax(i_mc)]
%      y=[ymin(i_mc)-ld,y,ymin(i_mc)-ld]
%      b=0
% elseif horizontal && x(end)< x(1) %% Shoreline horizontal & land on the left side (up)
%      x=[xmax(i_mc),x,xmin(i_mc)]
%      y=[ymax(i_mc)+ld,y,ymax(i_mc)+ld]
%      b=1%%
% elseif ~horizontal &&  y(end)< y(1)   %% vertical Shoreline  & land on the left side 
%     x=[xmax(i_mc)+ld,x,xmax(i_mc)+ld]
%      y=[ymin(i_mc),y,ymax(i_mc)]
%      b=1%%
% elseif ~horizontal && y(end)< y(1)  %% vertical Shoreline  & land on the right side 
%      x=[xmin(i_mc)-ld,x,xmin(i_mc)-ld]
%      y=[ymax(i_mc),y,ymin(i_mc)]
%      b=0
% end

xb=[x(end-1),x(end),x(1),x(2)];
yb=[y(end-1),y(end),y(1),y(2)];
   
plot(xb,yb,'k','linewidth',3);
hold on 

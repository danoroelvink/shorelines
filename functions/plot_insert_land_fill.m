function [x,y,xb,yb,xmax,ymax,xmin,ymin,nmax] = plot_insert_land_fill(x,y,ld,it,xmax,ymax,xmin,ymin,i_mc,n_mc,nmax,usefillpoints)
% function [x,y,xb,yb,xmax,ymax,xmin,ymin,nmax] = plot_insert_land_fill(x,y,ld,it,xmax,ymax,xmin,ymin,i_mc,n_mc,nmax,usefillpoints)
% 
% This routine defines the land-polygon behind the opne shoreline elements, 
% which are used for the plots. An approach can be used with 2 points at
% the landward side of the polygon (usefillpoints==0) or with a specified
% number of points landward of the coastline element (usefillpoints>=1). 
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2020 IHE Delft & Deltares
%
%       Dano Roelvink
%       d.roelvink@un-ihe.org
%       Westvest 7
%       2611AX Delft
%
%       Bas Huisman
%       bas.huisman@deltares.nl
%       Boussinesqweg 1
%       2629HV Delft
%
%   This library is free software: you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation, either
%   version 2.1 of the License, or (at your option) any later version.
%
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses>
%   --------------------------------------------------------------------
  
    dist=((y(end)-y(1)).^2+(x(end)-x(1)).^2).^0.5;
    distmax=((x(end)-x(1)).^2+(y(end)-y(1)).^2).^0.5; %max( max(((y-y(1)).^2+(x-x(1)).^2).^0.5), max(((y-y(end)).^2+(x-x(end)).^2).^0.5));
    ds0=hypot(x(2)-x(1),y(2)-y(1));
    
    %% determine orientation of shoreline
    if usefillpoints==0 || length(x)<3
        if it==0  % this value should store
            xmin(i_mc)=min(x);
            ymin(i_mc)=min(y);
            xmax(i_mc)=max(x);
            ymax(i_mc)=max(y);
            nmax=n_mc;
        elseif i_mc>length(xmin)
            xmin(i_mc)=min(x);
            ymin(i_mc)=min(y);
            xmax(i_mc)=max(x);
            ymax(i_mc)=max(y);
            nmax=n_mc;
        end
        o=atan2d((y(end)-y(1)),(x(end)-x(1)));
        
        % close the polygon
        if dist>0.2*distmax & usefillpoints==0
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
            xb=[x(end),x(1)];
            yb=[y(end),y(1)];
        else
            xb=[x(end-1),x(end),x(1),x(2)];
            yb=[y(end-1),y(end),y(1),y(2)];
        end
        
    else
        % number of points generated landward of coastline to close it
        % usefillpoints=6;
    
        % create id's of points where plotting points landward of open coast are generated
        % with ido for the orientation and idx for the x and y coordinates
        ido=ceil([1:(length(x)-2)/max(usefillpoints*2-2,1):length(x)-1]);
        idx=ceil([1:(length(x)-1)/(length(ido)-1):length(x)]);
        ido4=ido([1,2:2:length(ido)-1,length(ido)]);
        idx4=idx([1,2:2:length(idx)-1,length(idx)]);
        idx4=[idx4(1),floor((idx4(2:end-2)+idx4(3:end-1))/2),idx4(end)];
        o=mod(atan2d(diff(x),diff(y))-90,360);
        o2=[];
        for oo=1:usefillpoints
            o2(oo)=mean(get_smoothdata(o([ido4(oo):ido4(oo+1)]),'angleavgmean',3));
        end
        x2=x(idx4);
        y2=y(idx4);        
        
        % close the polygon
        if dist>0.2*distmax % may be check that dist>2*ds0 ?                %&& n_mc==1
            xb=[x(end),fliplr(x2-sind(o2)*ld),x(1)];
            yb=[y(end),fliplr(y2-cosd(o2)*ld),y(1)];
            x=[x,fliplr(x2-sind(o2)*ld),x(1)];
            y=[y,fliplr(y2-cosd(o2)*ld),y(1)];
        else
            xb=[x(end-1),x(end),x(1),x(2)];
            yb=[y(end-1),y(end),y(1),y(2)];
        end
    end
    
    %plot(xb,yb,'k-','linewidth',2);
    %hold on;plot(x,y,'g.','linewidth',3);
    %hold on 
    
end

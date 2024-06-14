function [xcr,ycr,indi,indj,ui,uj]=find_crossing(xi,yi,xj,yj)
% function [xcr,ycr,indi,indj,ui,uj]=find_crossing(xi,yi,xj,yj)
% 
% finds all crossings of polygons
% xcr and ycr are the crossing points.
% indi and indj provide the index of the vertex of the polygons i and j.
% (so indi=10 corresponds with points i=10 to i=11)
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
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses
%   --------------------------------------------------------------------

    xcr=[];
    ycr=[];
    indi=[];
    indj=[];
    ui=[];
    uj=[];
    nn=0;    
    xcr0=nan(length(xi)-1,length(xj)-1);
    ycr0=nan(length(xi)-1,length(xj)-1);
    eps=1e-12;
    for i=1:length(xi)-1
        for j=1:length(xj)-1
            dx1=diff(xi(i:i+1));
            dy1=diff(yi(i:i+1));
            dx2=diff(xj(j:j+1));
            dy2=diff(yj(j:j+1));

            if dx1~=0 && dx2~=0
                rc1 = dy1/dx1;
                rc2 = dy2/dx2;
                y1r = yi(i)-xi(i)*rc1;
                y2r = yj(j)-xj(j)*rc2;
                if (rc1-rc2)==0
                    % parallelle lijnen
                    xcr0(i,j)=nan;
                    ycr0(i,j)=nan;
                end
                xcr0(i,j) = (y2r-y1r)./(rc1-rc2);
                ycr0(i,j) = rc1*xcr0(i,j)+y1r;
                %ycr0 = rc2*xcr+y2r;
            elseif dx1==0 && dx2~=0
                rc2 = dy2/dx2;
                y2r = yj(j)-xj(j)*rc2;
                xcr0(i,j) = xi(i);
                ycr0(i,j) = rc2*xcr0(i,j)+y2r;
            elseif dx1~=0 && dx2==0
                rc1 = dy1/dx1;
                y1r = yi(i)-xi(i)*rc1;
                xcr0(i,j) = xj(j);
                ycr0(i,j) = rc1*xcr0(i,j)+y1r;
            elseif dx1==0 && dx2==0
                % parallelle lijnen
                xcr0(i,j)=nan;
                ycr0(i,j)=nan;
            end
            % check if it is on line segment          
            if xcr0(i,j)<max(min(xi(i:i+1)),min(xj(j:j+1)))-eps || xcr0(i,j)>min(max(xi(i:i+1)),max(xj(j:j+1)))+eps
                xcr0(i,j)=nan;
                ycr0(i,j)=nan;
            elseif ycr0(i,j)<max(min(yi(i:i+1)),min(yj(j:j+1)))-eps || ycr0(i,j)>min(max(yi(i:i+1)),max(yj(j:j+1)))+eps
                xcr0(i,j)=nan;
                ycr0(i,j)=nan;
            elseif ~isnan(xcr0(i,j))
                nn=nn+1;
                xcr(nn)=xcr0(i,j);
                ycr(nn)=ycr0(i,j);
                indi(nn)=i;
                indj(nn)=j;
                ui(nn)=((xcr(nn)-xi(i)).*(xi(i+1)-xi(i))+(ycr(nn)-yi(i))*(yi(i+1)-yi(i))) ./ ((xi(i+1)-xi(i)).^2+(yi(i+1)-yi(i)).^2);
                uj(nn)=((xcr(nn)-xj(j)).*(xj(j+1)-xj(j))+(ycr(nn)-yj(j))*(yj(j+1)-yj(j))) ./ ((xj(j+1)-xj(j)).^2+(yj(j+1)-yj(j)).^2);
            end
        end
    end
    
    [~,idu]=unique(indi+ui);
    xcr=xcr(idu);
    ycr=ycr(idu);
    indi=indi(idu);
    indj=indj(idu);
    ui=ui(idu);
    uj=uj(idu);
    [~,idu]=unique(indj+uj);
    xcr=xcr(idu);
    ycr=ycr(idu);
    indi=indi(idu);
    indj=indj(idu);
    ui=ui(idu);
    uj=uj(idu);
    
    %figure;plot(xi,yi,'b.-');hold on;plot(xj,yj,'r.-');
    %plot(xcr,ycr,'k*');
end


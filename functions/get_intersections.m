function [xcr,ycr,indc,inds,indi,indj,ui,uj]=get_intersections(xi0,yi0,xj0,yj0)
% function [xcr,ycr,indc,inds,indi,indj,ui,uj]=get_intersections(xi0,yi0,xj0,yj0)
%
% or
%
% function [xcr,ycr,indc,inds,indi,indj,ui,uj]=get_intersections(xi0,yi0,disentangle)
% 
% finds all crossings of polygons
% xcr and ycr are the crossing points.
% indi and indj provide the index of the vertex of the polygons i and j.
% (so indi=10 corresponds with points i=10 to i=11)
% 
% finds all crossings of polygons
% xcr and ycr are the crossing points.
% indi and indj provide the index of the vertex of the polygons i and j.
% (so indi=10 corresponds with points i=10 to i=11)
%
% INPUT:
%     xi        : x-coordinates of polygon 1 [m]
%     yi        : y-coordinates of polygon 1 [m]
%     xj        : x-coordinates of polygon 2 [m]
%     yj        : y-coordinates of polygon 2 [m]
%     disentangle : alternative 3rd argument which is used to disentangle xi/yi polygon 
% OUTPUT: 
%     xcr       : x-coordinates of crossings [m]
%     ycr       : y-coordinates of crossings [m]
%     indc      : index on polygon 1 of the crossing (=index of last point + fraction of last vertex) 
%     inds      : index on polygon 2 of the crossing (=index of last point + fraction of last vertex) 
%     indi      : index of last point of polygon 1 before the crossing (truncated version of indc)
%     indj      : index of last point of polygon 2 before the crossing (truncated version of inds)
%     ui        : fraction of polygon 1 up till crossing (from indi onwards)
%     uj        : fraction of polygon 2 up till crossing (from indj onwards)
% 
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2020 Deltares & IHE-Delft
%
%       Bas Huisman
%       bas.huisman@deltares.nl
%       Boussinesqweg 1
%       2629HV Delft
%
%       Dano Roelvink
%       d.roelvink@un-ihe.org
%       Westvest 7
%       2611AX Delft
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

    eps=1e-5;
    removesameindex=0;
    disentangle=0;
    if nargin==2
        xj0=xi0;
        yj0=yi0;
        removesameindex=1;
    elseif nargin==3
        xj0=xi0;
        yj0=yi0;
        removesameindex=1;
        disentangle=1;
    end
    xi0=xi0(:)';
    yi0=yi0(:)';
    xj0=xj0(:)';
    yj0=yj0(:)';
    
    xcr=[];
    ycr=[];
    indi=[];
    indj=[];
    ui=[];
    uj=[];
    
    m0=length(xi0)-1;
    n0=length(xj0)-1;
    
    xi=repmat(xi0(1:m0)',[1,n0]);
    yi=repmat(yi0(1:m0)',[1,n0]);
    xj=repmat(xj0(1:n0),[m0,1]);
    yj=repmat(yj0(1:n0),[m0,1]);
    
    ximin=min(repmat(xi0(1:m0)',[1,n0]),repmat(xi0(2:m0+1)',[1,n0]));
    ximax=max(repmat(xi0(1:m0)',[1,n0]),repmat(xi0(2:m0+1)',[1,n0]));
    yimin=min(repmat(yi0(1:m0)',[1,n0]),repmat(yi0(2:m0+1)',[1,n0]));
    yimax=max(repmat(yi0(1:m0)',[1,n0]),repmat(yi0(2:m0+1)',[1,n0]));

    xjmin=min(repmat(xj0(1:n0),[m0,1]),repmat(xj0(2:n0+1),[m0,1]));
    xjmax=max(repmat(xj0(1:n0),[m0,1]),repmat(xj0(2:n0+1),[m0,1]));
    yjmin=min(repmat(yj0(1:n0),[m0,1]),repmat(yj0(2:n0+1),[m0,1]));
    yjmax=max(repmat(yj0(1:n0),[m0,1]),repmat(yj0(2:n0+1),[m0,1]));
    
    nn=0;    
    dx1=repmat(diff(xi0)',[1,n0]);
    dy1=repmat(diff(yi0)',[1,n0]);
    dx2=repmat(diff(xj0),[m0,1]);
    dy2=repmat(diff(yj0),[m0,1]);

    % compute the derivative/gradient
    rc1 = dy1./dx1;
    rc2 = dy2./dx2;
    y1r = yi-xi.*rc1;
    y2r = yj-xj.*rc2;
    
    % compute potential crossings xcr1,xcr2 and xcr3
    xcr1 = (y2r-y1r)./(rc1-rc2);
    ycr1 = rc1.*xcr1+y1r;
    xcr2 = xi;
    ycr2 = rc2.*xcr2+y2r;
    xcr3 = xj;
    ycr3 = rc1.*xcr3+y1r;
    
    % apply potential crossings for cases where they can exist (e.g. where dx1 | dy1 is not zero)
    xc=nan(m0,n0);
    xc(dx1==0 & dx2~=0)=xcr2(dx1==0 & dx2~=0);
    xc(dx1~=0 & dx2==0)=xcr3(dx1~=0 & dx2==0);
    xc(dx1~=0 & dx2~=0)=xcr1(dx1~=0 & dx2~=0);
    yc=nan(m0,n0);
    yc(dx1==0 & dx2~=0)=ycr2(dx1==0 & dx2~=0);
    yc(dx1~=0 & dx2==0)=ycr3(dx1~=0 & dx2==0);
    yc(dx1~=0 & dx2~=0)=ycr1(dx1~=0 & dx2~=0);
    %cs(dx1==0 & dx2==0 | (rc1-rc2)==0)=4;

    % check if it is on line segment
    % remove any crossings that are beyond the length of the considered segments
    idnan1=xc<max(ximin,xjmin)-eps | xc>min(ximax,xjmax)+eps;
    xc(idnan1)=nan;
    yc(idnan1)=nan;
    idnan2=yc<max(yimin,yjmin)-eps | yc>min(yimax,yjmax)+eps;
    xc(idnan2)=nan;
    yc(idnan2)=nan;
    
    % remove the same index crossings in case similar lines are used (around the diagonal)
    % so all segment crossings with at least 1 exactly the same segment are nanned out
    if removesameindex==1
        idremove=[1:m0+1:m0*n0];
        xc(idremove)=nan;
        yc(idremove)=nan;
        idremove=[2:m0+1:m0*n0];
        xc(idremove)=nan;
        yc(idremove)=nan;
        idremove=[m0+1:m0+1:m0*n0];
        xc(idremove)=nan;
        yc(idremove)=nan;
        
        % make sure to use only the crossing points for the first line and not (the same mirrored ones) for the second line.
        if disentangle==0
        mat0=meshgrid([1:m0],[1:n0]);
        idmat=mat0>mat0';
        xc(idmat)=nan;
        yc(idmat)=nan;
        end
        
        % do not take beginning and end segment touching for a closed polygon
        if xi0(1)==xi0(end) && yi0(1)==yi0(end) && ~isempty(xc)
            xc(1,end)=nan;
            yc(1,end)=nan;
            xc(end,1)=nan;
            yc(end,1)=nan;
        end
    end
    
    % find the x,y coordinates of the crossing points 
    idnotnan=find(~isnan(xc));
    xcr=min(length(idnotnan),1);
    
    if nargout>1
        xcr=xc(idnotnan); 
        ycr=yc(idnotnan); 

        % find the indices of the line segment points just before the crossings (so never after!)
        [indi,indj]=find(~isnan(xc));
        
        % find the fraction of the line segement where the crossing is located 
        % (e.g. ui=0.25 means at 25% the length of the considered segment of xi,yi)
        % the indi and ui may be added to get the fraction of the polygon where the crossing is located 
        % (e.g. indj+uj = 16.3 means after the 16th node of xj,yj and at 30% of length in the direction of the 17th node of xj,yj)
        for nn=1:length(xcr)
            i=indi(nn);
            j=indj(nn);
            ui(1,nn)=((xcr(nn)-xi0(i)).*(xi0(i+1)-xi0(i))+(ycr(nn)-yi0(i))*(yi0(i+1)-yi0(i))) ./ ((xi0(i+1)-xi0(i)).^2+(yi0(i+1)-yi0(i)).^2);
            uj(1,nn)=((xcr(nn)-xj0(j)).*(xj0(j+1)-xj0(j))+(ycr(nn)-yj0(j))*(yj0(j+1)-yj0(j))) ./ ((xj0(j+1)-xj0(j)).^2+(yj0(j+1)-yj0(j)).^2);
        end
        
        % select only the unique points
        % and cosmetic changes to makes sure the vectors are horizontal
        ui=max(ui,0);
        uj=max(uj,0);
        ui=min(ui,1);
        uj=min(uj,1);
        [~,idu]=uniquetol(indi(:)+ui(:),eps);
        xcr=xcr(idu);
        ycr=ycr(idu);
        indi=indi(idu);
        indj=indj(idu);
        ui=ui(idu);
        uj=uj(idu);
        [~,idv]=uniquetol(indj(:)+uj(:),eps);
        xcr=xcr(idv);
        ycr=ycr(idv);
        indi=indi(idv);
        indj=indj(idv);
        ui=ui(idv);
        uj=uj(idv);
        %figure;plot(xi,yi,'b.-');hold on;plot(xj,yj,'r.-');
        %plot(xcr,ycr,'k*');
        
        xcr=double(xcr(:)');
        ycr=double(ycr(:)');
        indi=double(indi(:)');
        indj=double(indj(:)'); 
        ui=double(ui(:)');
        uj=double(uj(:)');
        indc=double(indi+ui);
        inds=double(indj+uj);
    end
end

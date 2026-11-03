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
    
    % Segment bounding boxes (column vectors for polygon i, row vectors for j)
    ximin=min(xi0(1:m0),xi0(2:m0+1))';
    ximax=max(xi0(1:m0),xi0(2:m0+1))';
    yimin=min(yi0(1:m0),yi0(2:m0+1))';
    yimax=max(yi0(1:m0),yi0(2:m0+1))';
    xjmin=min(xj0(1:n0),xj0(2:n0+1));
    xjmax=max(xj0(1:n0),xj0(2:n0+1));
    yjmin=min(yj0(1:n0),yj0(2:n0+1));
    yjmax=max(yj0(1:n0),yj0(2:n0+1));
    
    % Only segments whose bounding box overlaps the bounding box of the other
    % polygon can produce a crossing (the tolerance is larger than the 2*eps
    % used in the crossing test below), so the m0 x n0 crossing matrices are
    % only evaluated for those segments. This gives the same crossings as
    % evaluating all segment pairs, but is much faster for short lines.
    ri=(1:max(m0,0))';
    cj=1:max(n0,0);
    if removesameindex==0 && m0>0 && n0>0
        tol=4*eps;
        % segments next to a NaN (section separator), and vertical segments
        % when the other polygon has NaN-segments, are always kept, so that
        % results are exactly identical to the full evaluation
        nani=isnan(xi0(1:m0)+xi0(2:m0+1)+yi0(1:m0)+yi0(2:m0+1))';
        nanj=isnan(xj0(1:n0)+xj0(2:n0+1)+yj0(1:n0)+yj0(2:n0+1));
        nani=nani | (any(nanj) & (xi0(2:m0+1)-xi0(1:m0))'==0);
        nanj=nanj | (any(nani) & (xj0(2:n0+1)-xj0(1:n0))==0);
        ri=find(nani | (ximin<=max(xjmax)+tol & ximax>=min(xjmin)-tol & yimin<=max(yjmax)+tol & yimax>=min(yjmin)-tol));
        if isempty(ri)
            cj=find(nanj & false);
        else
            cj=find(nanj | (xjmin<=max(ximax(ri))+tol & xjmax>=min(ximin(ri))-tol & yjmin<=max(yimax(ri))+tol & yjmax>=min(yimin(ri))-tol));
        end
        ri=ri(:);
        cj=cj(:)';
    end
    ximin=reshape(ximin(ri),[],1);ximax=reshape(ximax(ri),[],1);yimin=reshape(yimin(ri),[],1);yimax=reshape(yimax(ri),[],1);
    xjmin=reshape(xjmin(cj),1,[]);xjmax=reshape(xjmax(cj),1,[]);yjmin=reshape(yjmin(cj),1,[]);yjmax=reshape(yjmax(cj),1,[]);
    m1=length(ri);
    n1=length(cj);
    
    xi=reshape(xi0(ri),[m1,1]);
    yi=reshape(yi0(ri),[m1,1]);
    xj=reshape(xj0(cj),[1,n1]);
    yj=reshape(yj0(cj),[1,n1]);
    dx1=reshape(xi0(ri+1),[m1,1])-xi;
    dy1=reshape(yi0(ri+1),[m1,1])-yi;
    dx2=reshape(xj0(cj+1),[1,n1])-xj;
    dy2=reshape(yj0(cj+1),[1,n1])-yj;

    % crossing of the two (infinite) lines through each pair of segments
    % (implicit expansion: column vectors for polygon i, row vectors for j)
    rc1 = dy1./dx1;
    rc2 = dy2./dx2;
    y1r = yi-xi.*rc1;
    y2r = yj-xj.*rc2;
    
    xc = (y2r-y1r)./(rc1-rc2);
    yc = rc1.*xc+y1r;
    both = (dx1~=0) & (dx2~=0);
    xc(~both)=nan;
    yc(~both)=nan;
    
    % vertical segment on polygon i
    id2 = (dx1==0) & (dx2~=0);
    if any(id2(:))
        xcr2 = xi+zeros(1,n1);
        ycr2 = rc2.*xcr2+y2r;
        xc(id2)=xcr2(id2);
        yc(id2)=ycr2(id2);
    end
    % vertical segment on polygon j
    id3 = (dx1~=0) & (dx2==0);
    if any(id3(:))
        xcr3 = xj+zeros(m1,1);
        ycr3 = rc1.*xcr3+y1r;
        xc(id3)=xcr3(id3);
        yc(id3)=ycr3(id3);
    end
    xc=reshape(xc,[m1,n1]);
    yc=reshape(yc,[m1,n1]);

    % remove crossings that lie outside the segments
    idnan1=xc<max(ximin,xjmin)-eps | xc>min(ximax,xjmax)+eps;
    xc(idnan1)=nan;
    yc(idnan1)=nan;
    idnan2=yc<max(yimin,yjmin)-eps | yc>min(yimax,yjmax)+eps;
    xc(idnan2)=nan;
    yc(idnan2)=nan;
    
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
    if nargout<=1
        % only whether there is a crossing (0 or 1)
        xcr=min(length(idnotnan),1);
        return
    end
    xcr=reshape(xc(idnotnan),[],1);
    ycr=reshape(yc(idnotnan),[],1);

    % find the indices of the line segment points just before the crossings (so never after!)
    [indi,indj]=find(~isnan(xc));
    indi=ri(indi);
    indj=cj(indj);
    indi=indi(:);
    indj=indj(:);
    
    % fraction along segment i and j of each crossing
    ui=(((xcr-xi0(indi)').*(xi0(indi+1)'-xi0(indi)')+(ycr-yi0(indi)').*(yi0(indi+1)'-yi0(indi)')) ./ ((xi0(indi+1)'-xi0(indi)').^2+(yi0(indi+1)'-yi0(indi)').^2))';
    uj=(((xcr-xj0(indj)').*(xj0(indj+1)'-xj0(indj)')+(ycr-yj0(indj)').*(yj0(indj+1)'-yj0(indj)')) ./ ((xj0(indj+1)'-xj0(indj)').^2+(yj0(indj+1)'-yj0(indj)').^2))';
    
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


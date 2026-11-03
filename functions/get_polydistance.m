function [dmin,xcr,ycr,ujmin,dmax,xcm,ycm,ujmax] = get_polydistance(Xr,Yr,Xc,Yc,Lcrit)
%function [dmin,xcr,ycr,ujmin,dmax,xcm,ycm,ujmax] = get_polydistance(Xr,Yr,Xc,Yc,Lcrit)
%
% INPUT: 
%    Xr          : x-coordinate of reference line [m]
%    Yr          : y-coordinate of reference line [m]
%    Xc          : x-coordinate of coastline [m]
%    Yc          : y-coordinate of coastline [m]
%
% OUTPUT:
%    dmin        : distance from reference line point to the specified coastline [m] (negative is at right-side of the line)
%    xcr         : x-coordinate of coastline point as projected on the grid [m] (closest value to the grid)
%    ycr         : y-coordinate of coastline point as projected on the grid [m] (closest value to the grid)
%    ujmin       : alongshore index of the cross-ing point for each grid point of the reference line [i.e. index as a fraction of Xc/Yc grid]
%    dmax        : distance from reference line point to the specified coastline for the farthest crossing [m] (negative is at right-side of the line)(farthest crossing from the grid)
%    xcm         : x-coordinate of coastline point as projected on the grid for the farthest crossing [m] (farthest crossing from the grid)
%    ycm         : y-coordinate of coastline point as projected on the grid for the farthest crossing [m] (farthest crossing from the grid)
%    ujmax       : alongshore index of the cross-ing point for each grid point of the reference line for the farthest crossing [i.e. index as a fraction of Xc/Yc grid]
%
% EXAMPLE:
%    Xr=[20.5:185.5]';
%    Yr=1-Xr/185;
%    Xc=[1:200]';
%    Yc=2+sin(Xc.^0.5).^2;
%    [dmin,xcr,ycr]=get_segmentdistance(Xr, Yr, Xc, Yc);
%
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2022 Deltares
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

    % dx,dy values
    if nargin<5
        Lcrit=500;  % maximum cross-shore distance to reference line (at both sides)
    end
    Xr=Xr(:);
    Yr=Yr(:);
    dx=Xr(2:end)-Xr(1:end-1); %for line in-between each reference point
    dy=Yr(2:end)-Yr(1:end-1);
    dx2 = [dx(1);(dx(2:end)+dx(1:end-1))/2;dx(end)];
    dy2 = [dy(1);(dy(2:end)+dy(1:end-1))/2;dy(end)];
    xcr=nan(length(Xr),1);
    ycr=nan(length(Xr),1);
    dmin=nan(length(Xr),1);
    ujmin=nan(length(Xr),1);
    if nargout>4
    xcm=nan(length(Xr),1);
    ycm=nan(length(Xr),1);
    dmax=nan(length(Xr),1);
    ujmax=nan(length(Xr),1);
    end
    Lxy=(dx2.^2+dy2.^2).^0.5;
    Xrnormal = [Xr,Xr]+[-Lcrit./Lxy.*dy2,Lcrit./Lxy.*dy2];
    Yrnormal = [Yr,Yr]+[Lcrit./Lxy.*dx2,-Lcrit./Lxy.*dx2];
    
    % Intersect all cross-shore normals with the coastline at once. Only
    % (normal, coastline segment) pairs whose bounding boxes overlap can cross
    % (pairs with NaN-adjacent segments are always kept, as get_intersections
    % can report those), so the crossing computation is done for those pairs
    % only. This uses the same arithmetic as get_intersections, and picks per
    % normal the crossing closest to the seaward end of the normal (smallest
    % ui, first coastline segment on ties), as the former loop did.
    eps=1e-5;                         % same tolerance as in get_intersections
    tol=4*eps;                        % bounding-box tolerance (> 2*eps)
    Xc=Xc(:);
    Yc=Yc(:);
    nseg=length(Xc)-1;
    if nseg>=1
        % normals (A) and coastline segments (B) as column vectors
        xa=Xrnormal(:,1); ya=Yrnormal(:,1);
        dxa=Xrnormal(:,2)-xa; dya=Yrnormal(:,2)-ya;
        xamin=min(Xrnormal,[],2); xamax=max(Xrnormal,[],2);
        yamin=min(Yrnormal,[],2); yamax=max(Yrnormal,[],2);
        nana=isnan(sum(Xrnormal,2)+sum(Yrnormal,2));
        xb=Xc(1:nseg); yb=Yc(1:nseg);
        dxb=Xc(2:end)-xb; dyb=Yc(2:end)-yb;
        xbmin=min(xb,Xc(2:end)); xbmax=max(xb,Xc(2:end));
        ybmin=min(yb,Yc(2:end)); ybmax=max(yb,Yc(2:end));
        nanb=isnan(xb+Xc(2:end)+yb+Yc(2:end));
        ov=(xamin<=xbmax'+tol & xamax>=xbmin'-tol & yamin<=ybmax'+tol & yamax>=ybmin'-tol) ...
           | nana | nanb';
        [ia,ib]=find(ov);
        if ~isempty(ia)
            xi=xa(ia); yi=ya(ia); dx1=dxa(ia); dy1=dya(ia);
            xj=xb(ib); yj=yb(ib); dx2=dxb(ib); dy2=dyb(ib);
            rc1=dy1./dx1;
            rc2=dy2./dx2;
            y1r=yi-xi.*rc1;
            y2r=yj-xj.*rc2;
            xc=(y2r-y1r)./(rc1-rc2);
            yc=rc1.*xc+y1r;
            both=(dx1~=0) & (dx2~=0);
            xc(~both)=nan;
            yc(~both)=nan;
            id2=(dx1==0) & (dx2~=0);
            xc(id2)=xi(id2);
            yc(id2)=rc2(id2).*xi(id2)+y2r(id2);
            id3=(dx1~=0) & (dx2==0);
            xc(id3)=xj(id3);
            yc(id3)=rc1(id3).*xj(id3)+y1r(id3);
            idnan=xc<max(xamin(ia),xbmin(ib))-eps | xc>min(xamax(ia),xbmax(ib))+eps;
            xc(idnan)=nan;
            yc(idnan)=nan;
            idnan=yc<max(yamin(ia),ybmin(ib))-eps | yc>min(yamax(ia),ybmax(ib))+eps;
            xc(idnan)=nan;
            yc(idnan)=nan;
            ok=~isnan(xc);
            ia=ia(ok); ib=ib(ok); xc=xc(ok); yc=yc(ok);
            xi=xi(ok); yi=yi(ok); dx1=dx1(ok); dy1=dy1(ok);
            xj=xj(ok); yj=yj(ok); dx2=dx2(ok); dy2=dy2(ok);
            
            % fraction along the normal (0 = seaward end, 1 = landward end)
            ui=((xc-xi).*dx1+(yc-yi).*dy1)./(dx1.^2+dy1.^2);
            ui=min(max(ui,0),1);
            uj=((xc-xj).*dx2+(yc-yj).*dy2)./(dx2.^2+dy2.^2);
            uj=min(max(uj,0),1);
            
            if ~isempty(ia)
                [~,order]=sortrows([ia,ui,ib]);
                first=order([true;diff(ia(order))~=0]);
                rows=ia(first);
                xcr(rows)=xc(first);
                ycr(rows)=yc(first);
                dmin(rows)=(0.5-ui(first))*Lcrit*2;
                ujmin(rows)=uj(first);              % alongshore location index
                if nargout>4
                    [~,order]=sortrows([ia,-ui,ib]);
                    first=order([true;diff(ia(order))~=0]);
                    rows=ia(first);
                    xcm(rows)=xc(first);
                    ycm(rows)=yc(first);
                    dmax(rows)=(0.5-ui(first))*Lcrit*2;
                    ujmax(rows)=uj(first);
                end
            end
        end
    end
end 

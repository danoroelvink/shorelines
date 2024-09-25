function [dmin,xcr,ycr,dmax,xcm,ycm] = get_polydistance(Xr,Yr,Xc,Yc,Lcrit)
%function [dmin,xcr,ycr,dmax,xcm,ycm] = get_polydistance(Xr,Yr,Xc,Yc,Lcrit)
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
%    dmin        : distance from reference line point to the specified coastline [m] (negative is at right-side of the line)(farthest crossing from the grid)
%    xcr         : x-coordinate of coastline point as projected on the grid [m] (farthest crossing from the grid)
%    ycr         : y-coordinate of coastline point as projected on the grid [m] (farthest crossing from the grid)
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
    if nargout>3
    xcm=nan(length(Xr),1);
    ycm=nan(length(Xr),1);
    dmax=nan(length(Xr),1);
    end
    Lxy=(dx2.^2+dy2.^2).^0.5;
    Xrnormal = [Xr,Xr]+[-Lcrit./Lxy.*dy2,Lcrit./Lxy.*dy2];
    Yrnormal = [Yr,Yr]+[Lcrit./Lxy.*dx2,-Lcrit./Lxy.*dx2];
    
    for ii=1:length(Xr)
        %quick=1;
        [xcr0,ycr0,~,~,~,~,ui,~]=get_intersections(Xrnormal(ii,:),Yrnormal(ii,:),Xc,Yc);
        if ~isempty(xcr0)
            imin=find(ui==min(ui));
            xcr(ii)=xcr0(imin);
            ycr(ii)=ycr0(imin);
            dmin(ii)=(0.5-ui(imin))*Lcrit*2;
            if nargout>3
                imax=find(ui==max(ui));
                xcm(ii)=xcr0(imax);
                ycm(ii)=ycr0(imax);
                dmax(ii)=(0.5-ui(imax))*Lcrit*2;
            end
        end 
    end 
end 

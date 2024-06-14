function [icm1,icm2,icp1,icp2]=find_neighbours(x_mc,y_mc,ic);
%% Find neighbouring points in x_mc,y_mc given index i and taking
%% into account cyclic or not
%   Detailed explanation goes here
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

    section=get_sectionnr(x_mc);
    [x, y, n_mc, i1, i2] = get_one_polygon(x_mc, y_mc, section(ic));
    n=length(x)-1;
    i=ic-i1+1;
    cyclic=hypot(x(end)-x(1),y(end)-y(1))<10;
    if cyclic
        im1=mod2(i-1,n);
        im2=mod2(i-2,n);
        ip1=mod2(i+1,n);
        ip2=mod2(i+2,n);
    else
        im1=max(i-1,1);
        im2=max(i-2,1);
        ip1=min(i+1,n+1);
        ip2=min(i+2,n+1);
    end
    icm1=im1+i1-1;
    icm2=im2+i1-1;
    icp1=ip1+i1-1;
    icp2=ip2+i1-1;
end

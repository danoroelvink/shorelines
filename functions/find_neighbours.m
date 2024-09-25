function [icm1,icm2,icp1,icp2]=find_neighbours(x_mc,y_mc,i_mc)
% function [icm1,icm2,icp1,icp2]=find_neighbours(x_mc,y_mc,i_mc)
%
% Find neighbouring points in x_mc,y_mc given index i and taking
% into account cyclic or not.
%
% INPUT: 
%    x_mc        : x-coordinates of grid [m]
%    y_mc        : y-coordinates of grid [m]
%    i_mc        : index of coastal element
% 
% OUTPUT:
%    icm1        : index -1
%    icm2        : index -2
%    icp1        : index +1
%    icp2        : index +2
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

    section=get_sectionnr(x_mc);
    [x, y, n_mc, i1, i2] = get_one_polygon(x_mc, y_mc, section(i_mc));
    n=length(x)-1;
    i=i_mc-i1+1;
    cyclic=hypot(x(end)-x(1),y(end)-y(1))<10;
    if cyclic
        im1=get_mod(i-1,n);
        im2=get_mod(i-2,n);
        ip1=get_mod(i+1,n);
        ip2=get_mod(i+2,n);
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

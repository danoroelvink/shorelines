function [ds0i,ds0i_center]=get_gridsize(x,y,ds0)
% function [dsi,dsi_center]=get_gridsize(x,y,ds0)
% 
% OUTPUT:
%         x       : y-coordinates of considered coastal section [m]
%         y       : number of points of considered coastal section
%         ds0     : grid size. Either a single ds0 value, or a [nx3] with columns for [xg,yg,ds0]
%
% Example:
%         dsi     : grid size for elements in-between coastal grid cells [n-1x1]
%         dsic    : grid size for coastline grid points [nx1]
%         dsiq    : grid size for transport grid points [n+1x1] (if xq and yq are available)
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

    if isscalar(ds0)
        ds0i=repmat(ds0,[1,length(x)]);
    else
        xg=ds0(:,1);
        yg=ds0(:,2);
        ds1=ds0(:,3);
        ds0i=get_interpolation_on_grid('weighted_distance',x,y,xg,yg,ds1);
    end
    ds0i_center=0.5*(ds0i(1:end-1)+ds0i(2:end));
    
    % remove double nans
    idnan=find(isnan(ds0i_center));
    idnotnan=setdiff([1:length(ds0i_center)],idnan(diff(idnan)==1));
    ds0i_center=ds0i_center(idnotnan);
    
    % transpose if needed
    ds0i=ds0i(:)';
    ds0i_center=ds0i_center(:)';
end

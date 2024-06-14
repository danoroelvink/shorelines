function [ cyclic ] = get_cyclic(x,y,ds0)
% function [ s ] = get_cyclic(x,y)
%
% Find out if polyline is cyclic; 1 if hypot(x(end-x(1),y(end)-y(1)))<frac*ds0;
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
    
    if ~isscalar(ds0)
        ds0=min(ds0(:,3));
    end
    frac   = 0.2;
    cyclic = hypot(x(end)-x(1),y(end)-y(1)) < frac*ds0;
end


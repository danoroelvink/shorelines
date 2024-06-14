function [distance, i]=get_disttopolyline(x,y,x0,y0,range)
%% Find shortest distance from point x0,y0 to polyline x,y
%    Syntax:
%    [distance,i]=locate_in_curvi_grid(x,y,x0,y0,range)
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

% This tool is part of <a href="http://OpenEarth.nl">OpenEarthTools</a>.
% OpenEarthTools is an online collaboration to share and manage data and 
% programming tools in an open source, version controlled environment.
% Sign up to receive regular updates of this function, and to contribute 
% your own tools.

    %% compute squared distance of x0,y0 to all points x,y where
    %% abs(x-x0)<range and abs(y-y0)<range
    inrange=abs(x-x0)<range&abs(y-y0)<range;
    r=zeros(size(x));
    r(inrange)=(x(inrange)-x0).^2+(y(inrange)-y0).^2;
    r(~inrange)=1.e10;
    %% Compute closest grid cell;
    %%
    [value,i]=min(r);
    distance=sqrt(value);
end
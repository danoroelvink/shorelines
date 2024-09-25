function [row,col]=find_in_grid(x,y,x0,y0);
% function [row,col]=find_in_grid(x,y,x0,y0);
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

    range=50000;
    inrange=abs(x-x0)<range&abs(y-y0)<range;
    dist=zeros(size(x));
    dist(inrange)=hypot(x(inrange)-x0,y(inrange)-y0);
    dist(~inrange)=1.e10;
    mindist=min(min(dist));
    [row,col,val]=find(dist==mindist);
end

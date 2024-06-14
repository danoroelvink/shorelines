function [coastS]=find_coastline(x,y,x_hard,y_hard)
% function [coastS]=find_coastline(x,y,x_hard,y_hard);
%
% UNTITLED Summary of this function goes here
% Detailed explanation goes here
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


    coastS=zeros(size(x_hard)-1);
    coastS=logical(coastS);
    for i=1:length(x_hard)-1
        [xx,yy]=get_intersections([x_hard(i),x_hard(i+1)],[y_hard(i),y_hard(i+1)],x,y);
        coastS(i)=int8(length(xx)>0);
    end
end

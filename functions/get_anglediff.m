function [angdif] = get_anglediff(angle1,angle2)
% function [angdif] = get_anglediff(angle1,angle2)
%
% Find the smallest angle (with sign) between two angles; all in degrees
% based on: cosd(angle)=u.v/(|u| |v|)
% with u=[cosd(angle1),sind(angle1)] and v=[cosd(angle2),sind(angle2)]
% Example:
%     get_anglediff(10,20)) returns -10
%     get_anglediff(10,350)) returns 20
% 
% angle1=round(angle1*1000)/1000;
% angle2=round(angle2*1000)/1000;
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
   
    % compute relative angle difference
    % choose smallest (-180° to +180°)
    angdif=mod(angle1-angle2+180,360)-180;
end


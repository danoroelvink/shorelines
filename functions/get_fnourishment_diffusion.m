function [K]=get_fnourishment_diffusion(V)
% function [K]=get_fnourishment_diffusion(V)
%
% Get the diffusion coefficient of a shoreface nourishment,
% which has been derived for the Dutch open coast based.
% 
% INPUT:
%    V    volume of the shoreface nourishment [m3]
%
% OUTPUT:
%    K    diffusion coefficient
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

    K = -4.86E-5 * V + 0.035 ; 
    K(K < 0.001) = 0.001; 

end 
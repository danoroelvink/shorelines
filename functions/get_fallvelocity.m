function [ws] = get_fallvelocity(D50, rhos, sal, temp)
%function [ws] = get_fallvelocity(D50, rhos, sal, temp)
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

    nu     = (1.14-0.031*(temp-15)+0.00068*(temp-15)^2)*10^-6;  % approximately 1e-6 at 20c;
    g      = 9.81;
    rhow   = max(1000+(0.8-(temp*0.002))*sal,1000);
    s      = rhos/rhow;
    Dster  = D50*((s-1)*g/(nu^2))^(1/3);
    ws     = nu/D50*(sqrt(10.36^2 + 1.049*Dster.^3)-10.36);

end

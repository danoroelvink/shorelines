function [kh,c] = wave_GUO2002(Tp,h)
% function [kh,c] = wave_GUO2002(Tp,h) 
%
% computes the kh number of the waves for shallow water
%
% INPUT:
%     Tp       peak wave period [s]
%     h        water depth [m]
%
% OUTPUT:
%     kh       wave number x water depth
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
    if Tp==0
	   error('GUO2002:TpGreaterThanZero','Tp must be greater than zero');
	end
    omega        = (2*pi)./Tp;
    kh           = omega.^2.*h/9.81.*(1-exp(-1*(omega.*(h/9.81).^0.5).^2.5)).^-0.4;
    c            = omega.*(h./kh);
end
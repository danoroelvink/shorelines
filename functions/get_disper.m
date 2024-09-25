function [k,C,Cg,n]=get_disper(h,T)
% function [k,C,Cg,n]=get_disper(h,T)
% 
% Approximate solution according to Guo (2002)
% 
% INPUT: 
%     h    : water depth [m]
%     T    : wave period [s]
% 
% OUTPUT:
%     k    : wave number [-]
%     C    : wave celerity [m/s]
%     Cg   : group velocity [m/s]
%     n    : deep/shallow water [0.5 to 1]
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2020 IHE Delft & Deltares
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

    g=9.81;
    sigma=2*pi./T;
    k = sigma.^2./g.*(1-exp(-(sigma.*sqrt(h./g)).^2.5)).^(-0.4);
    C = sigma./k;
    n = 0.5+k.*h./sinh(2.*k.*h);
    Cg=n.*C;
end

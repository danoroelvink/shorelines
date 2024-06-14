function plot_figureproperties(fh,HOR,VER,SCALE,X0,Y0)
% function plot_figureproperties(HOR,VER,SCALE,X0,Y0)
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

    if nargin<2;HOR=21*35;end
    if nargin<3;VER=29*35;end
    if nargin<4;SCALE=35;end
    if nargin<5;X0=70;end
    if nargin<6;Y0=70;end
    
    if (isoctave)
        pptyp = 'a4';
    else
        pptyp = 'a4letter';
    end    
    
    set(fh, 'PaperOrientation','portrait', ...
            'PaperPosition',[0 0 HOR/SCALE VER/SCALE], ...
            'PaperSize',[HOR/SCALE VER/SCALE], ...
            'Position',[X0,Y0,HOR,VER], ...
            'PaperUnits','centimeters', ...
            'PaperPositionMode','manual', ...
            'PaperType',pptyp); 
end

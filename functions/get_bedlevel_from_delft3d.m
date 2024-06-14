function [x,y,zb0,zb,dsdn]=get_bedlevel_from_delft3d(fname)
% function [x,y,zb0,zb,dsdn]=get_bedlevel_from_delft3d(fname)
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

    vs_use(fname);
    times=vs_let('map-info-series','ITMAPC');
    last=length(times);
    zb0=-squeeze(vs_let('map-sed-series',{1},'DPS','quiet'));
    zb=-squeeze(vs_let('map-sed-series',{last},'DPS','quiet'));
    dsdn=-squeeze(vs_let('map-const','GSQS','quiet'));
    sedero=zb-zb0;
    x=squeeze(vs_let('map-const','XZ','quiet'));
    y=squeeze(vs_let('map-const','YZ','quiet'));
    y(x==0)=nan;
    zb(x==0)=nan;
    x(x==0)=nan;
    if 1
        figure;
        pcolor(x,y,zb0)
        shading interp;
        axis equal;
        colorbar
    end
end    

function [MUD]=prepare_riverdischarges(S,COAST,TRANSP)
% function [MUD]=prepare_riverdischarges(S,COAST,TRANSP)
%
% This function prepares the river discharges and 
% initial conditions that are required to define mud supply by a river to the coast. 
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
    MUD=struct;
    if TRANSP.mud
        fprintf('  Prepare river discharges (mud) \n');
        MUD=struct;
        MUD.x=[];
        MUD.y=[];
        MUD.n=[];
        MUD.tstart=[];
        MUD.tend=[];
        MUD.rate_m3_per_yr=[];
    
        if ~isempty(findstr(lower(S.LDBriverdisch),'.riv'))
            % LBDriverdisch contains [xstart, ystart,  xend,  yend,  tstart(yyyymmdd), tend(yyyymmdd), totalvolume]
            riv=load(S.LDBriverdisch);
            MUD.n=size(riv,1);
            MUD.x=[riv(:,1),riv(:,3)];
            MUD.y=[riv(:,2),riv(:,4)];
            MUD.tstart=datenum(num2str(riv(:,5)),'yyyymmdd');
            MUD.tend=datenum(num2str(riv(:,6)),'yyyymmdd');
            MUD.rate_m3_per_yr=riv(:,7); % m3/year
        end
    end
end

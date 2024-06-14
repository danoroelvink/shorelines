function [MUD]=prepare_mudcoast(S)
% function [MUD]=prepare_mudcoast(S)
%
% This function prepares the use of the fine sediment properties and 
% initial conditions that are required to define a mud coast. 
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
    MUD.used   = S.mud;                                                                     % mud transport option
    if MUD.used
        % mud parameters
        MUD.taucr  = S.mud_taucr;
        MUD.M      = S.mud_M;
        MUD.B      = S.mud_B;
        MUD.w      = S.mud_w;
        MUD.Bfcrit = S.mud_Bfcrit;
        MUD.Bmmin  = S.mud_Bmmin;
        MUD.Bmmax  = S.mud_Bmmax;
        MUD.MHW    = S.mud_MHW;     % todo: check with climate change
        MUD.MSL    = S.mud_MSL;     % todo: check with climate change
        MUD.dm     = MUD.MHW-MUD.MSL;
        MUD.Tfm    = S.mud_Tfm;
        % river discharge initialisation
        MUD.xriv=[];
        MUD.yriv=[];
        MUD.nriv=[];
        MUD.tstart=[];
        MUD.tend=[];
        MUD.rate_m3_per_yr=[];
        % mangrove initialization
        MUD.xmgv=[];
        MUD.ymgv=[];
        MUD.nmgv=[];
        MUD.Bf=[];
        MUD.Bm=[];
        MUD.Bfm=[];
         
        if ~isempty(findstr(lower(S.LDBriverdisch),'.riv'))
            fprintf('  Prepare river discharges (mud) \n');
            % LBDriverdisch contains [xstart, ystart,  xend,  yend,  tstart(yyyymmdd), tend(yyyymmdd), volume/yr]
            riv=load(S.LDBriverdisch);
            MUD.nriv=size(riv,1);
            MUD.xriv=[riv(:,1),riv(:,3)];
            MUD.yriv=[riv(:,2),riv(:,4)];
            MUD.tstart=datenum(num2str(riv(:,5)),'yyyymmdd');
            MUD.tend=  datenum(num2str(riv(:,6)),'yyyymmdd');
            MUD.rate_m3_per_yr=riv(:,7); % m3/year
        end
        if ~isempty(findstr(lower(S.LDBmangrove),'.mgv'))
            fprintf('  Prepare mud flats and mangroves (mud) \n');
            % LBDmangrove contains [x,y,  Bf0,  Bm0]
            mgv=load(S.LDBmangrove);
            MUD.nmgv = size(mgv,1);
            MUD.xmgv = mgv(:,1);
            MUD.ymgv = mgv(:,2);
            MUD.Bf   = mgv(:,3); % mud flat width (m)
            MUD.Bm   = mgv(:,4); % mangrove width (m)
            MUD.Bfm  = mgv(:,5); % colonizing mangrove width (m)
        end
        %% write logfile
        % struct2log(MUD,'MUD','a');

    end
end

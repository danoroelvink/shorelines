function [MUD]=prepare_mudcoast(S)
% function [MUD]=prepare_mudcoast(S)
%
% This function prepares the use of the fine sediment properties and 
% initial conditions that are required to define a mud coast. 
%
% INPUT: 
%     S
%        .mud             : switch for using mud transport [0/1]
%        .mudtaucr        : critical shear stress for erosion [N/m2]
%        .mudm            : erosion rate [kg/m2/s]
%        .mudb            : muddy transport zone width [m]
%        .mudw            : fall velocity of the muddy sediment [m/s] 
%        .mudbfcrit       : critical mudflat width [m]
%        .mudbmmin        : minimum for the width of the muddy transport zone [m]
%        .mudbmmax        : maximum for the width of the muddy transport zone [m]
%        .mudmhw          : mhw level [m]
%        .mudmsl          : msl level [m]
%        .mudtfm          : tfm value
%        .ldbriverdisch   : river discharge file with .riv extension, containing [xriv1, yriv1, xriv2, yriv2, tstart, tend, rate]
%        .ldbmangrove     : mangrove definition file with .mgv extension, containing [xmgv, ymgv, Bf, Bm, Bfm] 
% 
% OUTPUT:
%     MUD
%        .used            : use mud transport [0/1]
%        .taucr           : critical shear stress for erosion [N/m2]
%        .M               : erosion rate [kg/m2/s]
%        .w               : fall velocity of the muddy sediment [m/s] 
%        .Bfcrit          : critical mudflat width [m]
%        .B               : muddy transport zone width [m]
%        .Bmmin           : minimum for the width of the muddy transport zone [m]
%        .Bmmax           : maximum for the width of the muddy transport zone [m]
%        .MHW             : mhw level [m]
%        .MSL             : msl level [m]
%        .dm              : high water level w.r.t. msl [m]
%        .Tfm             : tfm value
%        .ldbriverdisch   : river discharge file, containing [xriv1, yriv1, xriv2, yriv2, tstart, tend, rate]
%        .nriv            : number of the rivers
%        .xriv            : x-coordinates of rivers
%        .yriv            : y-coordinates of rivers
%        .tstart          : start time for mud supply [days in datenum format]
%        .tend		  : end time for mud supply [days in datenum format]
%        .rate		  : supply of mud by rivers [m3/yr]
%        .ldbmangrove     : mangrove definition file, containing [xmgv, ymgv, Bf, Bm, Bfm]
%        .nmgv            : number of the mangrove coasts
%        .xmgv            : x-coordinates of the mangrove coasts
%        .ymgv            : y-coordinates of the mangrove coasts
%        .Bf              : mud flat width [m]
%        .Bm  		  : mangrove width [m]
%        .Bfm 		  : colonizing mangrove width [m]
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
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses>
%   --------------------------------------------------------------------

    MUD=struct;
    MUD.used   = S.mud;                                                                     % mud transport option
    if MUD.used
        % mud parameters
        MUD.taucr  = S.mudtaucr;
        MUD.M      = S.mudm;
        MUD.B      = S.mudb;
        MUD.w      = S.mudw;
        MUD.Bfcrit = S.mudbfcrit;
        MUD.Bmmin  = S.mudbmmin;
        MUD.Bmmax  = S.mudbmmax;
        MUD.MHW    = S.mudmhw;
        MUD.MSL    = S.mudmsl;
        MUD.dm     = MUD.MHW-MUD.MSL;
        MUD.Tfm    = S.mudtfm;
        % river discharge initialisation
        MUD.xriv=[];
        MUD.yriv=[];
        MUD.nriv=[];
        MUD.tstart=[];
        MUD.tend=[];
        MUD.rate=[];
        % mangrove initialization
        MUD.xmgv=[];
        MUD.ymgv=[];
        MUD.nmgv=[];
        MUD.Bf=[];
        MUD.Bm=[];
        MUD.Bfm=[];
         
        if ~isempty(findstr(lower(S.ldbriverdisch),'.riv'))
            fprintf('  Prepare river discharges (mud) \n');
            % LBDriverdisch contains [xriv1, yriv1, xriv2, yriv2, tstart, tend, rate]
            riv=load(S.ldbriverdisch);
            MUD.nriv=size(riv,1);
            MUD.xriv=[riv(:,1),riv(:,3)];
            MUD.yriv=[riv(:,2),riv(:,4)];
            MUD.tstart=datenum(num2str(riv(:,5)),'yyyymmdd');
            MUD.tend=datenum(num2str(riv(:,6)),'yyyymmdd');
            MUD.rate=riv(:,7); % m3/year
        end
        if ~isempty(findstr(lower(S.ldbmangrove),'.mgv'))
            fprintf('  Prepare mud flats and mangroves (mud) \n');
            % LBDmangrove contains [x, y, Bf0,  Bm0, Bfm0]
            mgv=load(S.ldbmangrove);
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

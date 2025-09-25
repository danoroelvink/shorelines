function [WIND]=prepare_windconditions(S,TIME)
% function [WIND]=prepare_windconditions(S,TIME)
%
% This function reads and prepares the wind data.  
% A data-structure WIND is created, which is used throughout the computation. 
% 
% INPUT:
%   S   Structure with input data
%        .dune         : switch for computing dune evolution (0/1)
%        .mud          : switch for computing mud coast evolution (0/1)
%        .dt           : timestep of the dune computations [yr]
%        .rhoa         : density of air [kg/m3]
%        .Cd           : wind coefficient
%        .uz           : static value for the wind velocity (if not replaced by WND)
%        .phiwnd       : static value for the wind direction (if not replaced by WND)
%        .wndfile      : file of the wind input timeseries or climate, with extension .WND
%        .Windclimfile : alternative, filename (backup option)
%   TIME
%        .it           : number of timesteps since model start (it=0 at t0)
%
% OUTPUT: 
%   WIND   Structure with wind data
%        .dune         : switch for computing dune evolution (0/1)
%        .mud          : switch for computing mud coast evolution (0/1)
%        .dt           : timestep of the dune computations [yr]
%        .rhoa         : density of air [kg/m3]
%        .Cd           : wind coefficient
%        .uz           : static value for the wind velocity (if not replaced by WND)
%        .phiwnd       : static value for the wind direction (if not replaced by WND)
%        .wndfile      : file of the wind input timeseries or climate, with extension .WND
%        .WND
%           .uz        : timeseries of wind velocity [m/s] from wvcfile
%           .Dir       : timeseries of wind direction [°N] from wvcfile
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

    WIND         = struct;
    WIND.dune    = S.dune;
    WIND.mud     = S.mud;
    WIND.dt      = S.dtdune;
    WIND.mergeconditions=S.mergeconditions;

    WND=[];
    if S.dune || S.mud
        WIND.rhoa    = S.rhoa;
        WIND.Cd      = S.cd;
        WIND.uz      = S.uz;
        WIND.phiwnd  = S.phiwnd0;
        WIND.wndfile = S.wndfile;
        if isfield(S,'Windclimfile') && isempty(WIND.wndfile)
            WIND.wndfile = S.Windclimfile;
        end
        
        if ~isempty(WIND.wndfile)
            if ischar(WIND.wndfile)
                WIND.wndfile={WIND.wndfile};
            elseif size(WIND.wndfile,2)>3 && (size(WIND.wndfile,1)==1 || size(WIND.wndfile,1)==3)
                WIND.wndfile=WIND.wndfile';
            end
            
            % read wind data files
            fprintf('  Read wind conditions for runup at dunes\n');
            [WND]=get_inputfiledata(WIND.wndfile,TIME);
        end
        %% write logfile
        % struct2log(WIND,'WIND','a');
    end
    WIND.WND=WND;

end 

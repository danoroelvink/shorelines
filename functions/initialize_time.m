function [TIME]=initialize_time(S)
% function [TIME]=initialize_time(S)
% 
% INPUT:
%   S
%
% OUTPUT:
%   TIME
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2021 IHE Delft & Deltares
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
    fprintf('  Initialize time \n');

    TIME=struct;
    TIME.dt=S.dt;                                                              % timestep [years]
    TIME.tc=S.tc;                                                              % switch for using adaptive time step
    TIME.reftime=S.reftime;                                                    % Reference time (i.e. 'yyyy-mm-dd') <- leave empty to use t=0
    TIME.endofsimulation=S.endofsimulation;
    
    %% Reference start time of the model
    if ~isfield(S,'timenum0')
        if ~isempty(S.reftime)
            %timenum0 = datenum(S.reftime,'yyyy-mm-dd'); %HH:MM:SS
            timenum0 = datenum(S.reftime); %HH:MM:SS
        else
            timenum0 = 0;
        end
    else
        timenum0=S.timenum0;
    end
    S.times(1)=timenum0;
    TIME.timenum0=timenum0;
    
    %% Time parameters
    TIME.tnow=TIME.timenum0;
    if ~isfield(S,'tend')
        TIME.tend=datenum(S.endofsimulation);
    elseif isempty(S.tend)
        TIME.tend=datenum(S.endofsimulation); 
    else
        TIME.tend=S.tend;
    end
    TIME.tnext=TIME.tnow;  % next time for jpg output in 'make_Plot'
    TIME.it=-1;   % for now = -1 to follow the exisiting code
    TIME.itout=0;
    TIME.automatic=S.dt<=0; 
    TIME.tprev=TIME.tnow;

    %% write logfile
    % struct2log(TIME,'TIME','a');

end

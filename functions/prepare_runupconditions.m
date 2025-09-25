function [RUNUP]=prepare_runupconditions(S,TIME)
% function [RUNUP]=prepare_runupconditions(S,TIME)
%
% The runup data is initialized in this function. 
% Relevant input data is stored in the RUNUP data-structure. 
% 
% INPUT:
%    S
%         .swl0         : static value of the water-level [m w.r.t. MSL]
%         .watfile      : filename with time-series data of water-levels (leave empty to use static value)
%         .watclimfile  : alternative, older, filename for time-series data
%         .hso          : static value of the offshore wave height [m]
%         .tper		: static value of the offshore wave period [s]
%         .phiw0	: static value of the offshore wave direction [°N]
%         .wvdfile      : filename with time-series data of offshore waves (leave empty to use static values)
%         .waveclimfile : alternative, older, filename for time-series data
% 
% OUTPUT:
%    RUNUP
%         .dune         : switch for computing dune evolution (0/1)
%         .watfile      : filename with time-series data of water-levels (leave empty to use static value)
%         .swl0         : static value of the water-level [m w.r.t. MSL]
%         .SWL          : data structure with time-series of the water-level (with timenum and data field) from 'watfile', which overrules the static .swl
%             .timenum  : dates of water-level timeseries [days in datenum format]
%             .swl      : timeseries of surge water-levels [m]
%         .x            : field with the x-coordinates of the points with time-series of the water-level data
%         .y            : field with the y-coordinates of the points with time-series of the water-level data
%         .nloc         : number of points with time-series water-level data (is 0 when static values are used)
%         .wvdfile      : filename with time-series data of offshore waves (leave empty to use static values)
%         .Hs           : static value of the offshore wave height [m]
%         .Tp           : static value of the offshore wave period [s]
%         .Dir          : static value of the offshore wave direction [°N]
%         .WVD          : data structure with time-series of the offshore wave conditions (with timenum and data field) from 'wvdfile', which overrules the static .Hs, .Tp and .Dir
%             .timenum  : dates of wave timeseries [days in datenum format]
%             .Hs       : timeseries of offshore wave height [m]
%             .Tp       : timeseries of offshore wave period [s]
%             .Dir      : timeseries of offshore wave direction [°N]
%         .xw           : field with the x-coordinates of the points with time-series of the offshore wave data
%         .yw           : field with the y-coordinates of the points with time-series of the offshore wave data
%         .nlocw        : number of points with time-series offshore wave data (is 0 when static values are used)
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

    fprintf('  Prepare runup conditions \n');
    
    if S.dune || S.transmission==1
        
        %% Waterlevel information for dune evolution and wave transmission over breakwaters
        % create structure
        RUNUP=struct;        
        RUNUP.swl0=S.swl0;
        RUNUP.watfile=S.watfile;
        RUNUP.mergeconditions=S.mergeconditions;

        % Use regular wave climate of wvdfile is not specified
        if isempty(RUNUP.watfile)
            % backward compatibility with other input
            if isfield(S,'watclimfile')
            RUNUP.watfile=S.watclimfile;
            end
        end
        
        % Still water level data used for runup computations
        RUNUP.nloc=0;
        RUNUP.SWL=[];
        if ~isempty(RUNUP.watfile)   
            % Read input files for water levels
            fprintf(' Read waterlevel conditions for runup at dunes or wave transmission over breakwaters\n');
            [SWL]=get_inputfiledata(RUNUP.watfile,TIME);
            RUNUP.nloc=length(SWL);
            RUNUP.SWL=SWL;
            
            % Set the xy locations for the water level locations
            if isfield(S,'watlocfile')
                xyWat=load(S.watlocfile);
                RUNUP.x=xyWat(:,1);
                RUNUP.y=xyWat(:,2);
            else
                for kk=1:length(SWL)
                RUNUP.x(kk)=SWL(kk).x;
                RUNUP.y(kk)=SWL(kk).y;
                end
            end
        end

        if S.dune  

        %% Offshore wave information for runup 
            % create structure          
            RUNUP.dune=S.dune;
            RUNUP.Hs=S.hso;
            RUNUP.Tp=S.tper;
            RUNUP.Dir=S.phiw0;
            RUNUP.wvdfile=S.wvdfile;
    
            % Use regular wave climate of wvdfile is not specified
            if isempty(RUNUP.wvdfile)
                RUNUP.wvdfile=S.wvcfile;
                % backward compatibility with other input
                if isfield(S,'waveclimfile')
                    if ~isempty(S.waveclimfile)
                    RUNUP.wvdfile=S.waveclimfile;
                    end
                end
            end
            
            % Offshore wave data used for runup computations
            RUNUP.nlocw=0;
            RUNUP.WVD=[]; 
            if ~isempty(RUNUP.wvdfile)
                % Read input files for offshore waves
                fprintf('  Read offshore wave conditions for runup at dunes\n');
                [WVD]=get_inputfiledata(RUNUP.wvdfile,TIME);
                RUNUP.nlocw=length(WVD);
                RUNUP.WVD=WVD;
    
                % Set the xy locations for the wave locations
                if isfield(S,'WaveLocfile')
                    xywave=load(S.WaveLocfile);
                    RUNUP.xw=xywave(:,1);
                    RUNUP.yw=xywave(:,2);
                else
                    for kk=1:length(WVD)
                    RUNUP.xw(kk)=WVD(kk).x;
                    RUNUP.yw(kk)=WVD(kk).y;
                    end
                end
            end
                      
            % Check if the waterlevel and waves files are aligned
            if length(RUNUP.WVD)~=length(RUNUP.SWL)
                if (length(RUNUP.WVD)>1 && length(RUNUP.SWL)~=1) || (length(RUNUP.WVD)~=1 && length(RUNUP.SWL)>1)
                    fprintf('Warning : The number of waterlevel locations (.WAT-files) do not match the number of offshore wave locations (.WVD-files) \n')
                % elseif length(RUNUP.WVD)==1 && length(RUNUP.SWL)>1
                %     RUNUP.WVD=repmat(RUNUP.WVD,[length(RUNUP.SWL),1]);
                %     RUNUP.xw=RUNUP.x;
                %     RUNUP.yw=RUNUP.y;
                % elseif length(RUNUP.SWL)==1 && length(RUNUP.WVD)>1
                %     RUNUP.SWL=repmat(RUNUP.SWL,[length(RUNUP.WVD),1]);
                %     RUNUP.x=RUNUP.xw;
                %     RUNUP.y=RUNUP.yw;
                end
            end
        end
        %% write logfile
        % struct2log(RUNUP,'RUNUP','a');
    else
        RUNUP=[];
    end
end

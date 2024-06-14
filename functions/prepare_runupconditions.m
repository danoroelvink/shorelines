function [RUNUP]=prepare_runupconditions(S,TIME)
% function [RUNUP]=prepare_runupconditions(S,TIME)
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

    fprintf('  Prepare runup conditions \n');
    
    if S.dune || S.transmission==1
        
        %% Waterlevel information for dune evolution and wave transmission over breakwaters
        % create structure
        RUNUP=struct;        
        RUNUP.swl=S.SWL0;
        RUNUP.WATfile=S.WATfile;

        % Use regular wave climate of WVDfile is not specified
        if isempty(RUNUP.WATfile)
            % backward compatibility with other input
            if isfield(S,'Watclimfile')
            RUNUP.WATfile=S.Watclimfile;
            end
        end
        
        % Still water level data used for runup computations
        RUNUP.nloc=0;
        RUNUP.SWL=[];
        if ~isempty(RUNUP.WATfile)   
            % Read input files for water levels
            fprintf(' Read waterlevel conditions for runup at dunes or wave transmission over breakwaters\n');
            [SWL]=get_inputfiledata(RUNUP.WATfile,TIME);
            RUNUP.nloc=length(SWL);
            RUNUP.SWL=SWL;
            
            % Set the xy locations for the water level locations
            if isfield(S,'WatLocfile')
                xyWat=load(S.WatLocfile);
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
            RUNUP.Hs=S.Hso;
            RUNUP.Tp=S.tper;
            RUNUP.Dir=S.phiw0;
            RUNUP.WVDfile=S.WVDfile;
    
            % Use regular wave climate of WVDfile is not specified
            if isempty(RUNUP.WVDfile)
                RUNUP.WVDfile=S.WVCfile;
                % backward compatibility with other input
                if isfield(S,'Waveclimfile')
                    if ~isempty(S.Waveclimfile)
                    RUNUP.WVDfile=S.Waveclimfile;
                    end
                end
            end
            
            % Offshore wave data used for runup computations
            RUNUP.nlocw=0;
            RUNUP.WVD=[]; 
            if ~isempty(RUNUP.WVDfile)
                % Read input files for offshore waves
                fprintf('  Read offshore wave conditions for runup at dunes\n');
                [WVD]=get_inputfiledata(RUNUP.WVDfile,TIME);
                RUNUP.nlocw=length(WVD);
                RUNUP.WVD=WVD;
    
                % Set the xy locations for the wave locations
                if isfield(S,'WaveLocfile')
                    xyWave=load(S.WaveLocfile);
                    RUNUP.xw=xyWave(:,1);
                    RUNUP.yw=xyWave(:,2);
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

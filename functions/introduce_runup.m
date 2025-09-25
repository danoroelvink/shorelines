function [RUNUP]=introduce_runup(RUNUP,TIME,COAST,DUNE,STRUC,WAVE)
% function [RUNUP]=introduce_runup(RUNUP,TIME,COAST,DUNE,STRUC,WAVE)
% 
% Interpolate the relevant conditions that determine the 
% water-levels (run-up) at the dunes and transmissive offshore breakwaters 
% along the coastline at every timestep for multiple water-level 
% and offshore wave data sources. If only one wave source is available, 
% or when only static conditions are prescribed then no interpolation 
% is needed. Not only alongshore interpolation takes place [1xN], but  
% also the conditions are interpolated for each dune timestep (dtdune)  
% over the period of the coastal timestep (dt) [Mx1].
% 
% INPUT: 
%   RUNUP              : Structure with water-level and wave data of the ShorelineS model
%                        Containing the sub-structures for the water-levels and waves:
%        .dune         : switch for computing dune evolution (0/1)
%        .x            : field with the x-coordinates of the points with time-series water-level data
%        .y            : field with the y-coordinates of the points with time-series water-level data
%        .nloc         : number of points with time-series water-level data (is 0 when static values are used)
%        .swl          : static value of the water-level [m w.r.t. MSL]
%        .SWL          : data structure with time-series of the water-level (with timenum and data field) from 'watfile', which overrules .swl
%           .timenum   : dates of water-level timeseries [days in datenum format]
%           .swl       : timeseries of surge water-levels [m]
%        .xw           : field with the x-coordinates of the points with time-series offshore wave data
%        .yw           : field with the y-coordinates of the points with time-series offshore wave data
%        .nlocw        : number of points with time-series offshore wave data (is 0 when static values are used)
%        .Hs           : static value of the offshore wave height [m]
%        .Tp           : static value of the offshore wave period [s]
%        .Dir          : static value of the offshore wave direction [°N]
%        .WVD          : data structure with time-series of the offshore waves (with timenum and data field) from 'wvdfile', which overrules .swl
%           .Hs        : timeseries of offshore wave height [m]
%           .Tp        : timeseries of offshore wave period [s]
%           .Dir       : timeseries of offshore wave direction [°N]
%           .timenum   : dates of wave timeseries [days in datenum format]
%   TIME               : Structure with time/date information
%        .it           : timestep index (starts at it=0 at model start)
%        .tnow         : current time [in days in datenum format]
%        .tprev        : time at previous coastal timestep [in days in datenum format]
%   COAST              : Structure with x,y locations of coastline points
%   DUNE               : Dune information
%        .used         : switch for dune computations (0/1)
%        .dt           : timestep of the dune computations [yr]
%   STRUC              : Information of coastal structures
%        .transmission : switch for using wave transmission (0/1)
% 
% OUTPUT:
%   RUNUP              : Structure with water-level and wave data of the ShorelineS model + interpolated conditions
%        .swl          : interpolated value of the surge water-level [m w.r.t. MSL], [MxN] matrix with number of dune timesteps within next coastal timestep (M) times the coastline grid length (N)
%        .Hs           : interpolated value of the offshore wave height [m], [MxN] matrix with number of dune timesteps within next coastal timestep (M) times the coastline grid length (N)
%        .Tp           : interpolated value of the offshore wave period [s], [MxN] matrix with number of dune timesteps within next coastal timestep (M) times the coastline grid length (N)
%        .Dir          : interpolated value of the offshore wave direction [°N], [MxN] matrix with number of dune timesteps within next coastal timestep (M) times the coastline grid length (N) 
% 
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2024 IHE Delft & Deltares
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

    if DUNE.used || STRUC.transmission==1
        RUNUP.swl=RUNUP.swl0;
        if isempty(RUNUP)
            % static conditions specified as S.swl
            RUNUP.nloc=0;
        elseif isfield(RUNUP,'SWL')
            %% Interpolate conditions in time (all columns for all positions at once)
            x=COAST.x;
            y=COAST.y;
            
            %% interpolate SWL
            SWL=RUNUP.SWL;
            SWLnow=[];
            RUNUP.nloc=length(SWL);
            for kk=1:length(SWL)
                % interpolate swl at tnow
                eps=1e-6;          
                
                % check how many time instances of the SWL are within the current coastline timestep
                tnow=TIME.tnow+WAVE.dtrewind(1);
                tprev=TIME.tprev+WAVE.dtrewind(1);
                idt=[];
                if ~isempty(SWL(kk).timenum)
                idt=find(SWL(kk).timenum>tprev & SWL(kk).timenum<=tnow);
                end
                
                % construct dune timesteps in case a different dt is specified for the dunes than for the coastline
                timenow=tnow;
                if ~isempty(DUNE.dt) && TIME.it~=0
                    dt=min(DUNE.dt*365,tnow-tprev);
                    timenow=[tprev+dt:dt:tnow]';
                end

                if isfield(SWL,'htide') && ischar(RUNUP.swl0)
                    % using htide (second data column of the water-level file) in addition to the surge-data (first data column of the water-level file)
                    SWL(kk).htide=SWL(kk).htide;
                else
                    SWL(kk).htide=repmat(RUNUP.swl0,[length(SWL(1).swl),1]);
                end
                
                %% draw swl condition from an SWL-climate
                if ~isempty(WAVE.iwc) && isempty(SWL(kk).timenum)
                    % in case of a swl climate
                    if length(WAVE.WVC(1).Prob)==length(SWL(kk).Prob)
                        % use wave climate condition if SWL and WAVE are synchronized
                        RUNUP.iwc=WAVE.iwc;
                        RUNUP.Prob=WAVE.Prob;
                    elseif kk==1 && COAST.i_mc==1 && TIME.dtsteps==0
                        if RUNUP.mergeconditions==1
                            % draw all conditions with their respective probabilities
                            RUNUP.iwc=[1:length(SWL(1).swl)]';
                            RUNUP.Prob=SWL(1).Prob(:)./sum(SWL(1).Prob);
                        elseif ~isempty(SWL(1).Prob) 
                            % draw from swl climate based on probability
                            RUNUP.iwc=get_randsample(length(SWL(1).swl),1,SWL(1).Prob);
                            RUNUP.Prob=SWL(1).Prob(:)./sum(SWL(1).Prob);
                        else
                            % random number for drawing from swl climate
                            RUNUP.iwc=round((rand*length(SWL(1).swl)+0.5));
                            RUNUP.Prob=ones(length(SWL(1).swl),1)/length(SWL(1).swl);
                        end
                    end
                    SWLnow(:,kk) = SWL(kk).swl(RUNUP.iwc)+SWL(kk).htide(RUNUP.iwc);

                %% interpolate the water-levels from an SWL-timeseries
                elseif length(idt)>1 && min(abs(SWL(kk).timenum-tnow))<eps && min(abs(SWL(kk).timenum-tprev))<eps && isempty(DUNE.dt)
                    % in case more dense time-series than timesteps of coastline model
                    SWLnow(:,kk) = SWL(kk).swl(idt)+SWL(kk).htide(idt);
                elseif length(SWL(kk).timenum)>1 || TIME.it==0
                    % in case similar or larger timestep in data time-series than timesteps of coastline model
                    SWLnow(:,kk) = interp1(SWL(kk).timenum,SWL(kk).swl+SWL(kk).htide,timenow);
                else
                    % in case a single value is specified as the swl
                    SWLnow(1,kk) = SWL(kk).swl+SWL(kk).htide;
                end
            end
        end
        if RUNUP.nloc==0
            RUNUP.swl=RUNUP.swl(1)*ones(size(COAST.x));
        elseif RUNUP.nloc==1
            % Use uniform swl
            RUNUP.swl=SWLnow.*ones(size(SWLnow,1),length(COAST.x));
        else
            % Interpolate conditions along the coast
            method='weighted_distance';
            xw=RUNUP.x;
            yw=RUNUP.y;
            RUNUP.swl=get_interpolation_on_grid(method,x,y,xw,yw,SWLnow);
        end
    end
    if DUNE.used
        %% interpolate WVD
        WVD=RUNUP.WVD;
        RUNUP.nlocw=length(WVD);
        for kk=1:length(WVD)
            % interpolate swl at tnow
            eps=1e-6;          
            
            % check how many time instances of the SWL are within the current coastline timestep
            tnow=TIME.tnow+WAVE.dtrewind(1);
            tprev=TIME.tprev+WAVE.dtrewind(1);
            idt=[];
            if ~isempty(WVD(kk).timenum)
            idt=find(WVD(kk).timenum>tprev & WVD(kk).timenum<=tnow);
            end
            
            % construct dune timesteps in case a different dt is specified for the dunes than for the coastline
            timenow=tnow;
            if ~isempty(DUNE.dt) && TIME.it~=0
                dt=min(DUNE.dt*365,tnow-tprev);
                timenow=[tprev+dt:dt:tnow]';
            end
            
            %% draw swl condition from an SWL-climate
            if ~isempty(WAVE.iwc) && isempty(WVD(kk).timenum)
                % in case of a WVD climate
                if length(WAVE.WVC(1).Prob)==length(WVD(kk).Prob)
                    % use wave climate condition if WVD and WAVE are synchronized
                    RUNUP.iwc=WAVE.iwc;
                    RUNUP.Prob=WAVE.Prob;
                elseif kk==1 && COAST.i_mc==1 && TIME.dtsteps==0
                    if RUNUP.mergeconditions==1
                        % draw all conditions with their respective probabilities
                        RUNUP.iwc=[1:length(WVD(1).Hs)]';
                        RUNUP.Prob=WVD(1).Prob(:)./sum(WVD(1).Prob);
                    elseif ~isempty(WVD(1).Prob) 
                        % draw from WVD climate based on probability
                        RUNUP.iwc=get_randsample(length(WVD(1).Hs),1,WVD(1).Prob);
                        RUNUP.Prob=WVD(1).Prob(:)./sum(WVD(1).Prob);
                    else
                        % random number for drawing from WVD climate
                        RUNUP.iwc=round((rand*length(WVD(1).Hs)+0.5));
                        RUNUP.Prob=ones(length(WVD(1).Hs),1)/length(WVD(1).Hs);
                    end
                end
                HSnow(:,kk) = WVD(kk).Hs(RUNUP.iwc)';
                TPnow(:,kk) = WVD(kk).Tp(RUNUP.iwc)';
                DIRnow(:,kk) = WVD(kk).Dir(RUNUP.iwc)';

            %% interpolate the water-levels from an SWL-timeseries
            elseif length(idt)>1 && min(abs(WVD(kk).timenum-tnow))<eps && min(abs(WVD(kk).timenum-tprev))<eps && isempty(DUNE.dt)
                % in case more dense time-series than timesteps of coastline model
                HSnow(:,kk) = WVD(kk).Hs(idt);
                TPnow(:,kk) = WVD(kk).Tp(idt);
                DIRnow(:,kk)= WVD(kk).Dir(idt);
            elseif length(WVD(kk).timenum)>1 || TIME.it==0
                % in case similar or larger timestep in data time-series than timesteps of coastline model
                HSnow(:,kk) = interp1(WVD(kk).timenum,WVD(kk).Hs,timenow);
                TPnow(:,kk) = interp1(WVD(kk).timenum,WVD(kk).Tp,timenow);
                sinDIRnow   = interp1(WVD(kk).timenum,sind(WVD(kk).Dir),timenow);
                cosDIRnow   = interp1(WVD(kk).timenum,cosd(WVD(kk).Dir),timenow);
                DIRnow(:,kk)= mod(atan2d(sinDIRnow,cosDIRnow),360);
            else
                % in case a single value is specified as the WVD
                HSnow(1,kk) = WVD(kk).Hs;
                TPnow(1,kk) = WVD(kk).Tp;
                DIRnow(1,kk)= WVD(kk).Dir;
            end
        end
        if RUNUP.nlocw==0
            RUNUP.Hs=RUNUP.Hs(1)*ones(size(COAST.x));
            RUNUP.Tp=RUNUP.Tp(1)*ones(size(COAST.x));
            RUNUP.Dir=RUNUP.Dir(1)*ones(size(COAST.x));
        elseif RUNUP.nlocw==1
            RUNUP.Hs=HSnow(:)*ones(size(COAST.x));
            RUNUP.Tp=TPnow(:)*ones(size(COAST.x));
            RUNUP.Dir=DIRnow(:)*ones(size(COAST.x));
        else
            % Interpolate conditions along the coast
            method='weighted_distance';
            x=COAST.x;
            y=COAST.y;
            xw=RUNUP.xw;
            yw=RUNUP.yw;
            var1=[];
            var2=[];
            var1.Hs=HSnow;
            var1.Tp=TPnow;
            var2.Dir=DIRnow;
            [var1i,var2i]=get_interpolation_on_grid(method,x,y,xw,yw,var1,var2);
            RUNUP.Hs=var1i.Hs;
            RUNUP.Tp=var1i.Tp;
            RUNUP.Dir=var2i.Dir;
        end
    end
end
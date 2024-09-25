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
        if isempty(RUNUP)
            % static conditions specified as S.swl
            RUNUP.nloc=1;
        elseif ~isfield(RUNUP,'timenum')
            %% Interpolate conditions in time (all columns for all positions at once)
            x=COAST.x;
            y=COAST.y;
            
            %% interpolate SWL
            SWL=RUNUP.SWL;
            SWLnow=[];
            for kk=1:length(SWL)
                % interpolate swl at tnow
                eps=1e-6;          
                
                % check how many time instances of the SWL are within the current coastline timestep
                idt=find(SWL(kk).timenum>TIME.tprev & SWL(kk).timenum<=TIME.tnow);
                
                % construct dune timesteps in case a different dt is specified for the dunes than for the coastline
                timenow=TIME.tnow;
                if ~isempty(DUNE.dt) && TIME.it~=0
                    dt=min(DUNE.dt*365,TIME.tnow-TIME.tprev);
                    timenow=[TIME.tprev+dt:dt:TIME.tnow]';
                end
                
                % interpolate the water-levels
                if length(idt)>1 && min(abs(SWL(kk).timenum-TIME.tnow))<eps && min(abs(SWL(kk).timenum-TIME.tprev))<eps && isempty(DUNE.dt)
                    % in case more dense time-series than timesteps of coastline model
                    SWLnow(:,kk) = SWL(kk).swl(idt);
                elseif length(SWL(kk).timenum)>1 || TIME.it==0
                    % in case similar or larger timestep in data time-series than timesteps of coastline model
                    SWLnow(:,kk) = interp1(SWL(kk).timenum,SWL(kk).swl,timenow);
                else
                    % in case a single value is specified as the swl
                    SWLnow(1,kk) = SWL(kk).swl;
                end
            end
        elseif isfield(RUNUP,'swl')
            % get wave climate
            if length(WAVE.WVC(1).Hs)==length(SWL(kk).swl) && ~isempty(WAVE.iwc)
                % obtain random climate condition from wave climate 
                RUNUP.iwc=WAVE.iwc;
            elseif COAST.i_mc==1 && isfield(RUNUP,'Prob')
                % determine random climate condition based on probability
                RUNUP.iwc=get_randsample(length(RUNUP(1).swl),1,RUNUP(1).Prob);
            else
                % determine random climate condition assuming equal probability of conditions
                RUNUP.rnd=rand;
                RUNUP.iwc=round((RUNUP.rnd*length(RUNUP(1).swl)+0.5));
            end 
            % use a wind climate conditiosn 
            RUNUP.swl=[];
            for kk=1:length(RUNUP)
                RUNUP.swl=RUNUP(kk).swl(RUNUP.iwc);      
                RUNUP.swl=RUNUP.swl(:)';
                if isfield(RUNUP,'x')
                xw(1,kk)=RUNUP(kk).x;
                yw(1,kk)=RUNUP(kk).y;
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
        for kk=1:length(WVD)
            % interpolate swl at tnow
            eps=1e-6;          
            
            % check how many time instances of the SWL are within the current coastline timestep
            idt=find(WVD(kk).timenum>TIME.tprev & WVD(kk).timenum<=TIME.tnow);
            
            % construct dune timesteps in case a different dt is specified for the dunes than for the coastline
            timenow=TIME.tnow;
            if ~isempty(DUNE.dt) && TIME.it~=0
                dt=min(DUNE.dt*365,TIME.tnow-TIME.tprev);
                timenow=[TIME.tprev+dt:dt:TIME.tnow]';
            end
            
            % interpolate the wave conditions
            if length(idt)>1 && min(abs(WVD(kk).timenum-TIME.tnow))<eps && min(abs(WVD(kk).timenum-TIME.tprev))<eps && isempty(DUNE.dt)
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
            RUNUP.Hs=HSnow*ones(size(COAST.x));
            RUNUP.Tp=TPnow*ones(size(COAST.x));
            RUNUP.Dir=DIRnow*ones(size(COAST.x));
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
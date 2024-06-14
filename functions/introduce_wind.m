function [WIND]=introduce_wind(WIND,TIME,WAVE,COAST)
% function [WIND]=introduce_wind(WIND,TIME,WAVE,COAST)
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
    
    WND=WIND.WND;
    WVC=WAVE.WVC;
    xw=[];
    yw=[];
    if WIND.dune || WIND.mud
        if isempty(WND)
            % static conditions specified as S.uz and S.phiwnd0
            kk=1;
        elseif ~isempty(WND(1).timenum)
            % loop over wind time-series locations 
            WIND.uz=[];
            WIND.phiwnd=[];
            for kk=1:length(WND)
                eps=1e-6;           
                % check how many time instances of the SWL are within the current coastline timestep
                idt=find(WND(kk).timenum>TIME.tprev & WND(kk).timenum<=TIME.tnow);

                % construct dune timesteps in case a different dt is specified for the dunes than for the coastline
                timenow=TIME.tnow;
                if ~isempty(WIND.dt) && TIME.it~=0
                    dt=min(WIND.dt*365,TIME.tnow-TIME.tprev);
                    timenow=[TIME.tprev+dt:dt:TIME.tnow]';
                end
                
                % interpolate the wind conditions
                if length(idt)>1 && min(abs(WND(kk).timenum-TIME.tnow))<eps && min(abs(WND(kk).timenum-TIME.tprev))<eps && isempty(WIND.dt)
                    % in case more dense time-series than timesteps of coastline model
                    WIND.uz(:,kk) = WND(kk).uz(idt);
                    WIND.phiwnd(:,kk) = WND(kk).Dir(idt);     
                elseif length(WND(kk).timenum)>1 || TIME.it==0
                    % in case similar or larger timestep in data time-series than timesteps of coastline model
                    WIND.uz(:,kk) = interp1(WND(kk).timenum,WND(kk).uz,timenow);
                    WIND.phiwnd(:,kk) = mod(atan2d(interp1(WND(kk).timenum,sind(WND(kk).Dir),timenow),interp1(WND(kk).timenum,cosd(WND(kk).Dir),timenow)),360);     
                else
                    % in case a single value is specified as the WND
                    WIND.uz(1,kk) = WND(kk).uz;
                    WIND.phiwnd(1,kk) = WND(kk).Dir;
                end
                if isfield(WND,'x')
                xw(1,kk)=WND(kk).x;
                yw(1,kk)=WND(kk).y;
                end
            end        
            
        elseif isfield(WND,'uz')
            % get wave climate
            if length(WVC(1).Hs)==length(WND(1).uz)
                % obtain random climate condition from wave climate 
                WIND.iwc=WAVE.iwc;
            elseif COAST.i_mc==1 && isfield(WND,'Prob')
                % determine random climate condition based on probability
                WIND.iwc=get_randsample(length(WND(1).uz),1,WND(1).Prob);
            else
                % determine random climate condition assuming equal probability of conditions
                WIND.rnd=rand;
                WIND.iwc=round((WIND.rnd*length(WND(1).uz)+0.5));
            end 
            % use a wind climate conditiosn 
            WIND.uz=[];
            WIND.phiwnd=[];
            for kk=1:length(WND)
                WIND.uz=WND(kk).uz(WIND.iwc);      
                WIND.phiwnd=WND(kk).Dir(WIND.iwc);
                WIND.uz=WIND.uz(:)';
                WIND.phiwnd=WIND.phiwnd(:)';
                if isfield(WND,'x')
                xw(1,kk)=WND(kk).x;
                yw(1,kk)=WND(kk).y;
                end 
            end 
        end 
        if kk==1
            WIND.uz=WIND.uz(1)*ones(size(COAST.x));
            WIND.phiwnd=WIND.phiwnd(1)*ones(size(COAST.x));
        elseif kk>1
            % spatially varying wind
            xc=COAST.x;
            yc=COAST.y;
            method='weighted_distance';
            %method='alongshore_mapping';
            [WIND.uz,WIND.phiwnd]=get_interpolation_on_grid(method,xc,yc,xw,yw,WIND.uz,WIND.phiwnd);
        end         
    else
        WIND.phiwnd=[];
    end
end

%% do we need to include S.z ( elevation of measured data) here? or S.z should be manually inserted and used in dune evolution function?

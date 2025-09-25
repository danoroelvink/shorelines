function [TIME, WAVE] = get_coastlineupdate_timestep(TIME, BATHY, WAVE, FORMAT)
% function [TIME, WAVE] = get_coastlineupdate_timestep(TIME, BATHY, WAVE, FORMAT)
% 
% INPUT: 
%   TIME, BATHY, WAVE, FORMAT
%
% OUTPUT:
%   TIME, WAVE
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
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses>
%   --------------------------------------------------------------------

    %% Compute timestep of wave timeseries after current time
    WVCtime=1e9;
    if (~isempty(WAVE.wvcfile) || isfield(WAVE,'Wavematfile')) && ~isempty(WAVE.WVC(1).timenum)
        indxw=find((WAVE.WVC(1).timenum-TIME.tnow)>0,1);
        if ~isempty(indxw)
        WVCtime = WAVE.WVC(1).timenum(indxw);
        end
    end

    %% time step determination for coastline update
    eps=1e-5; % minimum time step
    tend=TIME.tend;
    if TIME.tnow+eps>tend
        tend=TIME.tnow+eps;
    end
    
    %% Cut timestep to output time point
    if (isempty(FORMAT.slplot) || FORMAT.tplot==0)
        FORMAT.tplot=1e9;
    end

    %% Cut timestep to bathy interaction time point
    if (isempty(BATHY.bathyupdate) || BATHY.tupdate==0) 
        BATHY.tupdate=1e9;
    end

    %% Cut timestep to flexible timestep ratio
    tc=1;
    if TIME.tc>0
        tc=TIME.tc;
    end

    %% Cut timestep to wave condition update time point / plot output point / bathy update point / end model time point
    dtnew=min([tc*TIME.adt,...
              (WVCtime-TIME.tnow)/365,...
              (tend-TIME.tnow)/365,...
              (FORMAT.tplot-TIME.tnow)/365,...
              (BATHY.tupdate-TIME.tnow)/365]);

    store.dt=TIME.dt;
    store.dtsteps=TIME.dtsteps;
    store.adt=TIME.adt;

    %% Set new timestep
    if TIME.tc>0
        % use dtnew when using time-varying timestep
        if isempty(WAVE.iwc)
            TIME.dt=dtnew;
            TIME.dtsteps=0;
        elseif TIME.dtsteps==0
            % repeat extreme wave conditions when a wave climate is used
            TIME.dt=dtnew;
            fac=ceil(TIME.dt0/(tc*TIME.adt));
            TIME.dtsteps=max(fac-1,0);
        else
            TIME.dtsteps=max(TIME.dtsteps-1,0);
        end

    elseif TIME.dtsteps>0 && TIME.tc==0
        % countdown the number of temprarily reduced timesteps when dtcrit<TIME.dt
        TIME.dtsteps=max(TIME.dtsteps-1,0);
        if TIME.dt>tc*TIME.adt
            % decrease timestep, but make sure it is an exact fraction of original dt
            fac=ceil(TIME.dt/(tc*TIME.adt));
            TIME.dt=TIME.dt/fac; 
            TIME.dtsteps=(TIME.dtsteps+1)*fac-1;
        elseif 2*TIME.dt<tc*TIME.adt && mod(TIME.dtsteps,2)==0
            % increase timestep, but make sure it exactly aligns with original dt
            dt2=(TIME.dtsteps+1)*TIME.dt; % remaining time
            fac=ceil(dt2/(tc*TIME.adt));
            TIME.dt=dt2/fac; 
            TIME.dtsteps=fac-1;
        end
    elseif TIME.dt0>(tc*TIME.adt) && TIME.tc==0 
        % use temporary smaller timestep ratio, 
        % but make sure it is an exact fraction of original dt
        fac=ceil(TIME.dt0/(tc*TIME.adt));
        TIME.dtsteps=fac-1;
        TIME.dt=TIME.dt0/fac; 
    else
        % use fixed timestep
        TIME.dt=TIME.dt0;
    end
    
    if 0
        % plot dt-steps
        if ~isfield(TIME,'dts');TIME.dts=[];end
        TIME.dts=[TIME.dts,TIME.dt];
        figure(11);clf;plot([0,cumsum(TIME.dts)],zeros(1,length(TIME.dts)+1),'b*');hold on;plot([0,cumsum(repmat(TIME.dts(1),[1 length(TIME.dts)]))],zeros(1,length(TIME.dts)+1),'ro');xlim([0,sum(TIME.dts)*1.01]);
    end
end


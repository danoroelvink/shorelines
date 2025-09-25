function [TIME]=get_nexttimestep(TIME,WAVE,DUNE,WIND)
% function [TIME]=get_nexttimestep(TIME,WAVE,DUNE,WIND)
% 
% This functions updates the timestep and prints the current 
% wave condition at each timestep of the simulation. 
%
% INPUT:
%     TIME      : TIME data structure with field tnow and trealstart
%     WAVE      : WAVE data structure with fields HStdp_mc and PHItdp_mc
%     DUNE      : DUNE data structure with field used
%     WIND      : WIND data structure with fields phiwnd, uz and Prob
% INPUT:
%     TIME      : TIME data structure with updated field tnow and tprev
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

    %% Next timestep
    TIME.tprev=TIME.tnow;
    TIME.tnow=TIME.tnow + TIME.dt*365;     %calculate the current time after each time step
    
    %% Display average time
    HSavg=median(WAVE.HStdp_mc(~isnan(WAVE.HStdp_mc)));
    PHIavg=mod(atan2d(median(sind(WAVE.PHItdp_mc(~isnan(WAVE.PHItdp_mc(:))))),median(cosd(WAVE.PHItdp_mc(~isnan(WAVE.PHItdp_mc(:)))))),360);
    if ~DUNE.used
        fprintf('   %s  Hs_tdp=%2.1f, Dir_tdp=%1.0f \n',datestr(double(TIME.tnow(1)),'yyyy-mm-dd HH:MM'),HSavg(1),PHIavg(1));     % fprintf('   %s \n',datestr(double(TIME.tnow),'yyyy-mm-dd HH:MM'));
    else
        wghtNR=1;
        if isfield(WIND,'Prob')
        wghtNR=repmat(WIND.Prob,[1,size(WIND.uz,2)]);
        end
        wghtWS=(WIND.uz.^1)/mean(sum(WIND.uz.^1,1));
        WSavg=sum(WIND.uz.*wghtWS.*wghtNR,1);
        WSavg=median(WSavg(~isnan(WSavg(:))));
        sdr = sum(sind(WIND.phiwnd).*wghtNR.*wghtWS,1);
        cdr = sum(cosd(WIND.phiwnd).*wghtNR.*wghtWS,1);
        WDavg=mod(atan2d(median(sdr(~isnan(sdr))),median(cdr(~isnan(cdr)))),360);
        fprintf('   %s  Hs_tdp=%2.1f, Dir_tdp=%1.0f, uz=%2.1f, uddir=%1.0f \n',datestr(double(TIME.tnow(1)),'yyyy-mm-dd HH:MM'),HSavg(1),PHIavg(1),WSavg(1),WDavg(1));  
    end
    
    %% Display model run-time
    if TIME.tnow>TIME.tend && ~(TIME.tnow==TIME.tend && TIME.tc==0)
        fprintf(['  Run completed successfully in ' num2str((now-TIME.trealstart)*86400) 's. \n']);
    end
    
end

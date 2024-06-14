function [CC]=prepare_climatechange(S, TIME)
% function [CC]=prepare_climatchange(S, TIME)
%
% Read climate change related corrections and initialize structure 
%
% INPUT:
%   S
%      .ccSLR    change in sea level, constant rate [m/yr] or time series of change since beginning simulation [m]
%      .ccHS     change in wave height constant rate [m/yr] or time series of change since beginning simulation [m]
%      .ccDIR    change in wave direction constant rate [deg/yr] or time series of change since beginning simulation [deg]
%
%
% OUTPUT:
%   CC
%       .timenum     time axis, datenum format
%       .SLR         instantaneous correction rate SLR (m)        
%       .HS          instantaneous correction rate Hs (m)
%       .DIR         instantaneous correction rate PHIw  (radians cartesian)
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2022 IHE Delft & Deltares
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

    fprintf('  Prepare climate change parameters \n');
    
    %% Initialize
    CC=struct;
    CC.timenum=NaN;    % serves as check whether we need introduce_climatechange or not
    CC.SLR=0.0;        % default rate of sea level rise [m/yr]
    CC.HS=0.0;         % default rate of increase in wave height [%/year]
    CC.DIR=0.0;        % default rate of wave direction change [°/yr]
    
    CC.HScor=1.0;      % default correction factor for the wave height [-], which is multiplied with the offshore wave height in 'introduce_wave' 
    CC.PHIcor=0.0;     % default correction factor for the wave direction [°], which is added to the offhore wave direction in 'introduce_wave'
    CC.SLRo=0.0;       % default rate of sea level rise [m/yr], which is used in the bruun-formula for coastal retreat in 'coastline_change' routine.
    
    %% 1. Sea level rise
    if isscalar(S.ccSLR)             
       CC.timenum=[];
       CC.SLR=S.ccSLR;   % constant rate ccSLR=0.002 [m/yr]
        
    elseif (ischar(S.ccSLR) && ~isempty(S.ccSLR)) || (isnumeric(S.ccSLR) && length(S.ccSLR)>=2)   % we have a timeseries with dMSL wrt initial MSL in time from file: yyyymmdd dMSL [m]
        if ischar(S.ccSLR)
            CCraw = load(S.ccSLR);
        else
            CCraw = S.ccSLR;
        end
        try
            CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmdd');
        catch
            try
                CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmddHHMM');
            catch
                CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmddHHMMSS');
            end
        end
        
        % Check length of input and abort if necessary
        if CC.timenum(1)>TIME.timenum0 || CC.timenum(end)<TIME.tend
           error('Time series for SLR not long enough'); 
        end
        
        CC.SLR=interpNANs(CCraw(:,2));      % correct derivative determined in introduce_climatechange 
        
    end
    
    
    %% 2. Wave characteristics correction
    if isscalar(S.ccHS)
        CC.timenum=[];
        CC.HS=S.ccHS; % constant rate of increase of HS [m/yr]
        
    elseif (ischar(S.ccHS) && ~isempty(S.ccHS)) || (isnumeric(S.ccHS) && length(S.ccHS)>=2)
        if ischar(S.ccHS)
            CCraw = load(S.ccHS);
        else
            CCraw = S.ccHS;
        end
        
        try
            CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmdd');
        catch
            try
                CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmddHHMM');
            catch
                CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmddHHMMSS');
            end
        end
        
        % Check length of input and abort if necessary
        if CC.timenum(1)>TIME.timenum0 || CC.timenum(end)<TIME.tend
           error('Time series for CC wave height change not long enough'); 
        end
        
        CC.HS=interpNANs(CCraw(:,2));    % relative change of HS with respect to initial situation at t0 [fraction]
        
    end
    
    if isscalar(S.ccDIR)
        CC.timenum=[];
        CC.DIR=S.ccDIR; % constant rate of change of the DIR [°/yr]
        
    elseif (ischar(S.ccDIR) && ~isempty(S.ccDIR)) || (isnumeric(S.ccDIR) && length(S.ccDIR)>=2)
        if ischar(S.ccDIR)
            CCraw = load(S.ccDIR);
        else
            CCraw = S.ccDIR;
        end
        try
            CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmdd');
        catch
            try
                CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmddHHMM');
            catch
                CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmddHHMMSS');
            end
        end
        
        % Check length of input and abort if necessary
        if CC.timenum(1)>TIME.timenum0 || CC.timenum(end)<TIME.tend
           error('Time series for CC wave direction change not long enough'); 
        end
        
        CC.DIR=interpNANsDIR(CCraw(:,2));    % relative change of DIR with respect to initial situation at t0 [fraction]
        
    end    
    %% write logfile
    % struct2log(CC,'CC','a');

end
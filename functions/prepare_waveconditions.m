function [WAVE]=prepare_waveconditions(S,TIME)
% function [WAVE]=prepare_waveconditions(S,TIME)
%
% This function reads and prepares the wave data.  
% A data-structure WAVE is created, which is used throughout the computation. 
% 
% OUTPUT:
%    WAVE
%      .ddeep                 : depth at the offshore boundary location for the waves (deep water) [m]
%      .dnearshore            : depth in the nearshore location [m] (e.g. at depth-of-closure)
%      .alpha                 : alpha parameter of the waves
%      .gamma                 : depth-induced wave breaking parameter
%      .waveinteraction       : wave interaction parameter
%      .tideinteraction       : tide interaction parameter
%      .surfwidth             : width of the surfzone [m]
%      .spread                : wave spreading [deg] (wave_dir from range:  S.phiw0 +/- 0.5*S.spread)
%      .diffraction           : switch for using diffraction (0/1)
%      .sphimax               : setting for the computation of the maximum angle of the waves (e.g. 42 degress org 'auto' which means that it is computed only at t0)
%      .interpolationmethod   : way of interpolating the wave conditions over the grid 
%      .wavefile              : file specifying the wave time-series and climate (specifies a wavegrid and gridded wave height/period, less common method) 
%      .wvcfile               : file specifying the wave time-series and climates (measured at observation points)
%      .WVC                   : data-structure with wave data that is read from the WVC files
%      .spacevaryingwave      : switch for space-varying wave conditions
%      .Hso                   : offshore wave height [m]
%      .phiw0                 : direction of the waves [°N]
%      .dirspr                : directional spreading paramter [°] at a moment in time in the 2d wave spectrum
%      .spread                : the spreading of the waves over time [°] which is not the directional spreading in the 2d wave spectrum
%      .tper                  : wave period [s]
%      .ccslr                 : climate impacted rise in sea level (SLR) [Nx2] with 'time in datenum format' and 'sea level with respect to initial situation' (rates per year are computed automatically) % S.tanbeta is used as 'slope angle'.
%      .ccdir                 : climate impacted change in wave direction [Nx2] with 'time in datenum format' and 'relative change in wave direction w.r.t. initial situation' as a # degrees (rates per year are computed automatically) 
%      .cchs                  : climate impacted change in wave height [Nx2] with 'time in datenum format' and 'relative change in wave height w.r.t. initial situation' as a non-dimensionless multiplicationfactor (rates per year are computed automatically) 
%      .DA                    : data-assimilation information
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
   
    fprintf('  Prepare wave conditions \n');
    WAVE=struct;
    WAVE.ddeep=S.ddeep;
    WAVE.dnearshore=S.dnearshore;
    WAVE.alpha=S.alpha;
    WAVE.gamma=S.gamma;
    WAVE.waveinteraction=S.waveinteraction;
    WAVE.tideinteraction=S.tideinteraction;
    WAVE.wavefile=S.wavefile;
    WAVE.surfwidth=S.surfwidth;
    WAVE.spread=S.spread;                                                             % wave spreading [?] (wave_dir from range:  S.phiw0 +/- 0.5*S.spread)
    WAVE.wvcfile=S.wvcfile;
    WAVE.diffraction=S.diffraction;    
    WAVE.diffsmooth=S.diffsmooth;
    WAVE.sphimax=S.sphimax;
    WAVE.interpolationmethod=S.interpolationmethod;  
    WAVE.mergeconditions=S.mergeconditions;

    % add climate impact
    WAVE.ccslr=S.ccslr;               % climate impacted rise in sea level (SLR) [Nx2] with 'time in datenum format' and 'sea level with respect to initial situation' (rates per year are computed automatically) % S.tanbeta is used as 'slope angle'.
    WAVE.ccdir=S.ccdir;               % climate impacted change in wave direction [Nx2] with 'time in datenum format' and 'relative change in wave direction w.r.t. initial situation' as a # degrees (rates per year are computed automatically) 
    WAVE.cchs=S.cchs;                 % climate impacted change in wave height [Nx2] with 'time in datenum format' and 'relative change in wave height w.r.t. initial situation' as a non-dimensionless multiplicationfactor (rates per year are computed automatically) 
    try WAVE.DA=S.da;end
    
    % add directional spreading
    if ~isfield(S,'dirspr')
        WAVE.dirspr=12;
    elseif isempty(S.dirspr)
        WAVE.dirspr=12;
    else
        WAVE.dirspr=S.dirspr;
    end
    
    % read wave input files, backward compatibility with old keywords
    WAVE.iwc=[];
    if isempty(WAVE.wvcfile) && isfield(S,'waveclimfile') 
        WAVE.wvcfile=S.waveclimfile;
    end
    
    % read content of files
    WAVE.WVC=[];
    if ~isempty(WAVE.wvcfile)
        % Read data from wave input files
        % various data-files can be read
        [WVC]=get_inputfiledata(WAVE.wvcfile,TIME);
        
        % Aggregate wave time-series data
        % method of speeding up the simulation by computing average conditions for time step periods
        [WVC]=get_aggregatedtimeseries(WVC,TIME);
        
        % Add wave time-series
        WAVE.WVC=WVC;
        WAVE.spacevaryingwave=WVC(1).spacevaryingwave;
    else
        WAVE.Hso=S.hso;
        WAVE.phiw0=S.phiw0;
        WAVE.spread=S.spread;
        WAVE.tper=S.tper;               
    end
    WAVE.dtrewind=zeros(1,max(length(WAVE.WVC),1));

end

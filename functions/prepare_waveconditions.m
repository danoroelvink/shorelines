function [WAVE]=prepare_waveconditions(S,TIME)
% function [WAVE]=prepare_waveconditions(S,TIME)
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
   
    fprintf('  Prepare wave conditions \n');
    WAVE=struct;
    WAVE.ddeep=S.ddeep;
    WAVE.dnearshore=S.dnearshore;
    WAVE.alpha=S.alpha;
    WAVE.gamma=S.gamma;
    WAVE.wave_interaction=S.wave_interaction;
    WAVE.tide_interaction=S.tide_interaction;
    WAVE.wavefile=S.wavefile;
    WAVE.surf_width=S.surf_width;
    WAVE.surf_width_w=S.surf_width_w;
    WAVE.spread=S.spread;                                                             % wave spreading [?] (wave_dir from range:  S.phiw0 +/- 0.5*S.spread)
    WAVE.WVCfile=S.WVCfile;
    WAVE.diffraction=S.diffraction;    
    WAVE.sphimax=S.sphimax;
    WAVE.interpolationmethod=S.interpolationmethod;  
    
    % add climate impact
    WAVE.ccSLR=S.ccSLR;               % climate impacted rise in sea level (SLR) [Nx2] with 'time in datenum format' and 'sea level with respect to initial situation' (rates per year are computed automatically) % S.tanbeta is used as 'slope angle'.
    WAVE.ccDIR=S.ccDIR;               % climate impacted change in wave direction [Nx2] with 'time in datenum format' and 'relative change in wave direction w.r.t. initial situation' as a # degrees (rates per year are computed automatically) 
    WAVE.ccHS=S.ccHS;                 % climate impacted change in wave height [Nx2] with 'time in datenum format' and 'relative change in wave height w.r.t. initial situation' as a non-dimensionless multiplicationfactor (rates per year are computed automatically) 
    try WAVE.DA=S.DA;end
    
    % add directional spreading
    if ~isfield(S,'dirspr')
        WAVE.dirspr=12;
    elseif isempty(S.dirspr)
        WAVE.dirspr=12;
    else
        WAVE.dirspr=S.dirspr;
    end
    
    % read wave input files, backward compatibility with old keywords
    if isempty(WAVE.WVCfile) && isfield(S,'Waveclimfile') 
        WAVE.WVCfile=S.Waveclimfile;
    end
    
    % read content of files
    WAVE.WVC=[];
    if ~isempty(WAVE.WVCfile)
        % Read data from wave input files
        % various data-files can be read
        [WVC]=get_inputfiledata(WAVE.WVCfile,TIME);
        
        % Aggregate wave time-series data
        % method of speeding up the simulation by computing average conditions for time step periods
        [WVC]=get_aggregatedtimeseries(WVC,TIME);
        
        % Add wave time-series
        WAVE.WVC=WVC;
        WAVE.spacevaryingwave=WVC(1).spacevaryingwave;
    else
        WAVE.Hso=S.Hso;
        WAVE.phiw0=S.phiw0;
        WAVE.spread=S.spread;
        WAVE.tper=S.tper;               
    end
    %% write logfile
    % struct2log(WAVE,'WAVE','a');

end

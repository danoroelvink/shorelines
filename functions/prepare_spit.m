function [SPIT]=prepare_spit(S)
% [SPIT]=prepare_spit(S)
%
% INPUT:
%    S
%      .spit_method                                         % overwash formulation
%      .spit_width                                          % width of tip of spit (used for overwash)
%      .spit_headwidth                                      % width of tip of spit (used for upwind correction)
%      .OWscale                                             % scales the rate of the overwash per timestep (i.e. what part of the deficit is moved to the backbarrier)
%      .Dsf                                                 % underwater part of active height for shoreface -> used only in spit-width function
%      .Dbb                                                 % underwater part of active height for back-barrier -> used only in spit-width function
%      .Bheight                                             % berm height used for overwash funciton (i.e. added to Dsf or Dbb)
%      .tide_interaction
%      .wave_interaction
%      .wavefile                                            % wave table (.mat)
%      .surf_width_w                                        % width of surf zone, where to collect the wave conditions from wave table
%      .surf_width                                          % width of surf zone, where to update the bathymetry
%      .bathy_update                                        % the dates when the bathymetry should be updated, the input should be in dates form, can accept more than one  {'yyyy-mm-dd'}
%         
% OUTPUT:
%    SPIT
%      .method                                              % overwash formulation (currently 'default')
%      .spit_width                                          % width of tip of spit (used for overwash)
%      .spit_headwidth                                      % width of tip of spit (used for upwind correction)
%      .OWscale                                             % scales the rate of the overwash per timestep (i.e. what part of the deficit is moved to the backbarrier)
%      .Dsf                                                 % underwater part of active height for shoreface -> used only in spit-width function
%      .Dbb                                                 % underwater part of active height for back-barrier -> used only in spit-width function
%      .Bheight                                             % berm height used for overwash funciton (i.e. added to Dsf or Dbb)
%      .tide_interaction
%      .wave_interaction
%      .wavefile                                            % wave table (.mat)
%      .surf_width_w                                        % width of surf zone, where to collect the wave conditions from wave table
%      .surf_width                                          % width of surf zone, where to update the bathymetry
%      .bathy_update                                        % the dates when the bathymetry should be updated, the input should be in dates form, can accept more than one  {'yyyy-mm-dd'}
%
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

    fprintf('  Prepare spit \n');
    %% Migrating spit
    SPIT=struct;
    SPIT.used=1;   
    SPIT.method=S.spit_method;                                                             % overwash formulation
    SPIT.spit_width=S.spit_width;                                                         % width of tip of spit (used for overwash)
    SPIT.spit_headwidth=S.spit_headwidth;                                                 % width of tip of spit (used for upwind correction)
    SPIT.OWscale=S.OWscale;                                                               % scales the rate of the overwash per timestep (i.e. what part of the deficit is moved to the backbarrier)
    SPIT.OWtimescale=S.OWtimescale;                                                       % timescale for the overwash process in years(i.e. what part of the deficit is moved to the backbarrier)
    % SPIT.OWtimescale overrides SPIT.OWscale if it is >0
    if SPIT.OWtimescale>0
        SPIT.OWscale=0;
    end
    SPIT.Dsf=S.spit_Dsf;                                                                     % underwater part of active height for shoreface -> used only in spit-width function
    SPIT.Dbb=S.spit_Dbb;                                                                  % underwater part of active height for back-barrier -> used only in spit-width function
    SPIT.Bheight=S.Bheight;                                                               % berm height used for overwash funciton (i.e. added to Dsf or Dbb)
    SPIT.tide_interaction=S.tide_interaction;
    SPIT.wave_interaction=S.wave_interaction;
    SPIT.wavefile=S.wavefile;                                                             % wave table (.mat)
    SPIT.surf_width_w=S.surf_width_w;                                                     % width of surf zone, where to collect the wave conditions from wave table
    SPIT.surf_width=S.surf_width;                                                         % width of surf zone, where to update the bathymetry
    SPIT.bathy_update=S.bathy_update;                                                     % the dates when the bathymetry should be updated, the input should be in dates form, can accept more than one  {'yyyy-mm-dd'};
    %% write logfile
    % struct2log(SPIT,'SPIT','a');

end

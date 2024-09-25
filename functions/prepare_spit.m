function [SPIT]=prepare_spit(S)
% function [SPIT]=prepare_spit(S)
%
% Initializes the information needed to compute spit development, and creates SPIT data-structure.
% 
% INPUT: 
%    S
%      .spitmethod         : overwash formulation
%      .spitwidth          : width of tip of spit (used for overwash)
%      .spitheadwidth      : width of tip of spit (used for upwind correction)
%      .owscale            : scales the rate of the overwash per timestep (i.e. what part of the deficit is moved to the backbarrier)
%      .Dsf                : underwater part of active height for shoreface -> used only in spit-width function
%      .Dbb                : underwater part of active height for back-barrier -> used only in spit-width function
%      .bheight            : berm height used for overwash funciton (i.e. added to Dsf or Dbb)
%      .tideinteraction    : tide interaction parameter
%      .waveinteraction    : wave interaction parameter
%      .wavefile           : wave table (.mat)
%      .surfwidthw         : width of surf zone, where to collect the wave conditions from wave table
%      .surfwidth          : width of surf zone, where to update the bathymetry
%      .bathyupdate        : the dates when the bathymetry should be updated, the input should be in dates form, can accept more than one  {'yyyy-mm-dd'}
%         
% OUTPUT:
%    SPIT
%      .method             : overwash formulation (currently 'default')
%      .spitwidth          : width of tip of spit (used for overwash)
%      .spitheadwidth      : width of tip of spit (used for upwind correction)
%      .owscale            : scales the rate of the overwash per timestep (i.e. what part of the deficit is moved to the backbarrier)
%      .Dsf                : underwater part of active height for shoreface -> used only in spit-width function
%      .Dbb                : underwater part of active height for back-barrier -> used only in spit-width function
%      .bheight            : berm height used for overwash funciton (i.e. added to Dsf or Dbb)
%      .tideinteraction    : tide interaction parameter
%      .waveinteraction    : wave interaction parameter
%      .wavefile           : wave table (.mat)
%      .surfwidth          : width of surf zone, where to update the bathymetry
%      .bathyupdate        : the dates when the bathymetry should be updated, the input should be in dates form, can accept more than one  {'yyyy-mm-dd'}
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
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses>
%   --------------------------------------------------------------------

    fprintf('  Prepare spit \n');
    %% Migrating spit
    SPIT=struct;
    SPIT.used=1;   
    SPIT.method=S.spitmethod;                                                         % overwash formulation
    SPIT.spitwidth=S.spitwidth;                                                       % width of tip of spit (used for overwash)
    SPIT.spitheadwidth=S.spitheadwidth;                                               % width of tip of spit (used for upwind correction)
    SPIT.owscale=S.owscale;                                                           % scales the rate of the overwash per timestep (i.e. what part of the deficit is moved to the backbarrier)
    SPIT.owtimescale=S.owtimescale;                                                   % timescale for the overwash process in years(i.e. what part of the deficit is moved to the backbarrier)
    % SPIT.owtimescale overrides SPIT.owscale if it is >0
    if SPIT.owtimescale>0
        SPIT.owscale=0;
    end
    SPIT.Dsf=S.spitdsf;                                                               % underwater part of active height for shoreface -> used only in spit-width function
    SPIT.Dbb=S.spitdbb;                                                               % underwater part of active height for back-barrier -> used only in spit-width function
    SPIT.bheight=S.bheight;                                                           % berm height used for overwash funciton (i.e. added to Dsf or Dbb)
    SPIT.tideinteraction=S.tideinteraction;
    SPIT.waveinteraction=S.waveinteraction;
    SPIT.wavefile=S.wavefile;                                                         % wave table (.mat)
    SPIT.surfwidth=S.surfwidth;                                                       % width of surf zone, where to update the bathymetry
    SPIT.bathyupdate=S.bathyupdate;                                                   % the dates when the bathymetry should be updated, the input should be in dates form, can accept more than one  {'yyyy-mm-dd'};
    %% write logfile
    % struct2log(SPIT,'SPIT','a');

end

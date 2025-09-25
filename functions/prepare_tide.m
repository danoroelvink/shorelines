function [TIDE]=prepare_tide(S)
% function [TIDE]=prepare_tide(S)
% 
% The tide input is initialized, and stored in the TIDE data-structure. 
% 
% INPUT: 
%    S                 : Structure with input information of the ShorelineS model, with relevant fields:
%         .tidefile    : File with alongshore distribution of tidal
% 
% OUTPUT:
%    TIDE
%        .x_stat       : x coordinates of tide support points (m)
%        .y_stat       : y coordinates of tide support points (m)
%        .eta_stat     : M2 and M4 water level amplitudes (m)
%        .detads_stat  : longshore gradients of M2 and M4 tidal amplitudes
%        .phi_stat     : phase (deg) of M2 and M4 water level components
%        .k_stat       : longshore wave number M2 and M4 (rad/m)
%        .ss_stat      : longshore mean surface slope driving residual current
%        .cf           : roughness factor [-]
%        .hmin         : minimum depth [m]
%        .hclosure     : depth-of-closure [m]
%        .x            : x-coordinates of cross-shore profiles [m]
%        .zb           : z-coordinates of cross-shore profiles [m]
%        .Ttide        : duration of a spring-neap cycle (745*60)
%        .nT           : number of points in a single tide (24)
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

    fprintf('  Prepare tide \n');
    if isfield(S,'tidefile') & ~isempty(S.tidefile)
        td                      = load(S.tidefile);
        TIDE.x_stat             = td(:,1);
        TIDE.y_stat             = td(:,2);
        TIDE.eta_stat(:,1:2)    = td(:,3:4);
        TIDE.detads_stat(:,1:2) = td(:,5:6);
        TIDE.phi_stat(:,1:2)    = td(:,7:8);
        TIDE.k_stat(:,1:2)      = td(:,9:10);
        TIDE.ss_stat            = td(:,11);
        TIDE.cf                 = S.cf;
        TIDE.hmin               = S.hmin;
        TIDE.hclosure           = S.hclosure;
        if isfield(S,'profile')
            S.tideprofile=S.profile;
        end
        xz                      = load(S.tideprofile);
        TIDE.x  = [xz(1,1):S.tidedx:xz(end,1)];
        TIDE.zb = interp1(xz(:,1),xz(:,2),TIDE.x);
        TIDE.Ttide              = 745*60;
        TIDE.nT                 = 24;
        TIDE.n                  = S.tiden; 
        %% write logfile
        % struct2log(TIDE,'TIDE','a');
    else
        TIDE=[];
    end
    
end

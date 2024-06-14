function [TIDE]=interpolate_tide(TIDE,COAST)
% function [TIDE]=interpolate_tide(TIDE,COAST)
%
% INPUT :
%      TIDE
%         .x_stat      : x coordinates of tide support points (m)
%         .y_stat      : y coordinates of tide support points (m)
%         .eta_stat    : M2 and M4 water level amplitudes (m)
%         .detads_stat : longshore gradients of M2 and M4 tidal amplitudes
%         .phi_stat    : phase (deg) of M2 and M4 water level components
%         .k_stat      : longshore wave number M2 and M4 (rad/m)
%         .ss_stat     : longshore mean surface slope driving residual
%                        current
%      COAST            
%         .xq          : x-coordinate of transport points (m)
%         .yq          : y-coordinate of transport points (m)
% 
% OUTPUT :
%      TIDE
%         .eta         : M2 and M4 water level amplitudes (m)
%         .detads      : longshore gradients of M2 and M4 tidal amplitudes
%         .phi         : phase (deg) of M2 and M4 water level components
%         .k           : longshore wave number M2 and M4 (rad/m)
%         .ss          : longshore mean surface slope driving residual
%                        current
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

    if ~isempty(TIDE)
        %fprintf('  Interpolate tide alongshore \n');
        var1.eta1    =   TIDE.eta_stat(:,1);
        var1.eta2    =   TIDE.eta_stat(:,2);
        var1.detads1 =   TIDE.detads_stat(:,1);
        var1.detads2 =   TIDE.detads_stat(:,2);
        var1.k1      =   TIDE.k_stat(:,1);
        var1.k2      =   TIDE.k_stat(:,2);
        var1.ss      =   TIDE.ss_stat;
        var2.phi1    =   TIDE.phi_stat(:,1);
        var2.phi2    =   TIDE.phi_stat(:,2);
        xtide        =   TIDE.x_stat;
        ytide        =   TIDE.y_stat;
        xq           =   COAST.xq;
        yq           =   COAST.yq;
        method       =   'weighted_distance';
        [var1i,var2i,idGRID]=get_interpolation_on_grid(method,xq,yq,xtide,ytide,var1,var2);
        eta(:,1)    = var1i.eta1;
        eta(:,2)    = var1i.eta2;
        detads(:,1) = var1i.detads1;
        detads(:,2) = var1i.detads2;
        k(:,1)      = var1i.k1;
        k(:,2)      = var1i.k2;
        ss          = var1i.ss;
        phi(:,1)    = var2i.phi1;
        phi(:,2)    = var2i.phi2;
        TIDE.eta    = eta;
        TIDE.detads = detads;
        TIDE.k      = k;
        TIDE.ss     = ss;
        TIDE.phi    = phi;
    end
    
end

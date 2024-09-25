function [FNOUR] = get_fnourishment(TIME,COAST,FNOUR,WAVE,i_mc)
% function [FNOUR] = get_fnourishment(TIME,COAST,FNOUR,WAVE,i_mc)
% 
% Computes the instantaneous source term of the shoreface nourishment to
% the coastline and the remaining volume of the shoreface nourishment
%
% INPUT:  
%    TIME               : time data structure, with fields
%    COAST              : coast data structure, with fields   
%    FNOUR              : shoreface nourishment data structure, with fields 
%    WAVE               : wave data structure, with fields
%    i_mc               : index of the active coastal element
%
% OUTPUT:
%     FNOUR
%       .q_grid_tot     : actual nourishment rate in [m3/m/yr]
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

eps=1e-6;  % default minimum volume for shoreface nourishment
FNOUR.q_tot = zeros(size(COAST.x)); % default initialization value of 0 for the rate at which shoreface nourishments nourish the coast

if FNOUR.fnourish == 1 % if shoreface nourishments are activated
    t = TIME.it+1;
        
    for n = 1:FNOUR.n % for each shoreface nourishment

        % Get relevant segments and fraction
        xn       = FNOUR.x(n,:);
        yn       = FNOUR.y(n,:);
        w        = FNOUR.w(n);
        mb       = FNOUR.mb;
        labda0   = FNOUR.labda0;
        K        = FNOUR.K(n);
        V0       = FNOUR.V(n);
        [fraction,i_mc_rel,idx,dxnour]   = get_fnourishment_fraction(xn,yn,COAST.x_mc,COAST.y_mc,COAST.n_mc,COAST.ds_mc,COAST.shadow_mc);
        FNOUR.i_mc_rel{n}(t,:) = i_mc_rel; % needed for joining the computed shoreface nourishment volumes of the coastal sections
        FNOUR.idx{t,i_mc,n}    = idx{i_mc};   % needed for plotting 
        FNOUR.q{t,i_mc}(n,:)   = zeros(size(COAST.x)); % initialize the rate of supply by the shoreface nourishment to the coast at 0
        
        if ismember(i_mc,i_mc_rel) % only for relevant segments, otherwise q stays zero as initiated above 

            if TIME.tnow >= FNOUR.tstart(n) && FNOUR.tstart(n) <= TIME.timenum0 && TIME.nt == 0
                % Initialize the volume of the shoreface nourishment that is present at the model starttime (t0)
                % in case the nourishment has started before t0
                V0dens    = V0 / dxnour(i_mc);
                idx_mc    = idx{i_mc};     
                dt        = 1 * 24 * 3600; % timestep used for the foreshore nourishment evolution [s]
                delta_t   = (TIME.timenum0 - FNOUR.tstart(n)) * 24 * 3600; % time between nourishment start and start of modeltime [in s]
                dtsteps   = ceil(delta_t/dt); % number of steps
                [~,~,Fac] = get_fnourishment_rate(dt, dtsteps, V0dens, labda0, w, mb, K, WAVE, idx_mc, xn, yn);

                FNOUR.Vt(n,t) = V0 * Fac;
                FNOUR.Vt_mc{i_mc}(n,t) = V0 * Fac * fraction(i_mc);
          
            elseif TIME.tnow >= FNOUR.tstart(n) && FNOUR.Vt(n,t) == 0
                % Initialize the volume of the shoreface nourishment when it occurs after the model starttime (t0)
                % Set starting volume at t1 or later (find solution for when tstart = timenum0) 
                FNOUR.Vt(n,t) = V0;
                FNOUR.Vt_mc{i_mc}(n,t) = V0 * fraction(i_mc);
                
            elseif TIME.tnow < FNOUR.tstart(n) 
                % set Vt_mc to 0 before start time of nourishment
                FNOUR.Vt(n,t) = 0;
                FNOUR.Vt_mc{i_mc}(n,t) = 0; 

            end
            
            % Active nourishment requiring computation of the volume decrease of Vt_mc (check if remaining volume > 0)
            if TIME.tnow > FNOUR.tstart(n) && FNOUR.Vt(n,t) > 10*eps
                % Nourishment erodes : Compute volume based on previous time step and supply to the coast
                V1          = FNOUR.Vt(n,t) * fraction(i_mc);
                V1dens      = V1 / dxnour(i_mc); % [m3/m]
                idx_mc      = idx{i_mc};
                dt          = TIME.dt * 365 * 24 * 60 * 60; % yrs -> [s]
                dtsteps     = 1;
                [V2dens,Qend] = get_fnourishment_rate(dt, dtsteps, V1dens, labda0, w, mb, K, WAVE, idx_mc, xn, yn);
                
                FNOUR.q{t,i_mc}(n,idx{i_mc}) = Qend * 60 * 60 * 24 * 365; % set the rate of supply by the shoreface nourishment to the coast [m3/m/yr]
                FNOUR.Vt_mc{i_mc}(n,t+1) = V2dens * dxnour(i_mc);         % the new volume of the shoreface nourishment after supply to the coast [m3]
                
                % Set the nourishment volume to a very small value ~= 0 if it is below a threshold
                % the model knows that the nourishment has already fully eroded if it finds this 'epsilon' volume. 
                FNOUR.Vt_mc{i_mc}(n,t+1) = max(FNOUR.Vt_mc{i_mc}(n,t+1),eps); 
                               
            else
                % No erosion : Use same volume as for previous time step
                FNOUR.Vt_mc{i_mc}(n,t+1) = FNOUR.Vt_mc{i_mc}(n,t);
            end

        end
        
        % Add all sources together
        FNOUR.q_tot = sum(FNOUR.q{t,i_mc},1);    
    end 
        
    %% Join shoreface nourishment volumes 
    % compute new volume of the elements by adding up remaining sediment for each section on which t has been applied
    if i_mc==COAST.n_mc
        for n = 1 : FNOUR.n % nr of nourishments 
            temp=zeros(1,max(FNOUR.i_mc_rel{n}(end,:)));
            for j = FNOUR.i_mc_rel{n}(end,:) % nr of relevant segments 
                temp(j)  =  FNOUR.Vt_mc{j}(n,t+1);
            end 
            FNOUR.Vt(n,t+1) = sum( temp ); 
        end 
    end
    
    %% store data in mc-fields
    if i_mc==1
        FNOUR.q_tot_mc=FNOUR.q_tot;
    else
        FNOUR.q_tot_mc=[FNOUR.q_tot_mc,nan,FNOUR.q_tot];
    end
    
end 


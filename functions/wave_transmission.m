function [Kt_total] = wave_transmission(x_hard,y_hard,COAST,SWL,Hs,Tp,BW_depth,BW_cheight,BW_width,BW_slope,BW_transm_form,BW_D50)

%% wave transmission calculation
%
% INPUT:
%     Hs             significant wave height [m]
%     Tp             peak wave period [s]
%     Dir            wave direction [deg]
%     SWL             water level [m]      
%     BW_depth       depth at breakwater location [m]
%     BW_cheight     breakwater crest height [m]
%     BW_width       breakwater crest width [m]
%     BW_slope       breakwater slope [-]
%     BW_transm_form used formulation for calculating transmission 
%     BW_D50         D50 of breakwater armour material [m]
%
% OUTPUT:
%     Kt           transmitted wave [0-1]
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2020 IHE Delft & Deltares
%
%       Dano Roelvink, Ahmed Elghandour
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

%% Initialize fluxes
    [ x_struc,y_struc,n_hard,~,~ ] = get_one_polygon(x_hard,y_hard,1);
    Kt = zeros(length(Hs),length(n_hard));
    Kt_total = zeros(length(Hs),length(n_hard));
    
    % Use sub-timesteps when multiple measurements of SWL are taking place in a single timestep 
    % (i.e. at more than just one timestep, but also inbetween moments)
    dtfractions=max([size(Hs,1),size(Tp,1),size(SWL,1)]);
  
    % Loop over the sub-conditions of the time step (i.e. at fractions of coastline timestep)
    for tt=1:dtfractions   
        
        %% Prepare SWL format
        SWLi=SWL(min(tt,size(SWL,1)),:);   
        % Interpolate SWL conditions along coast for wave locations
        method='weighted_distance';
        [SWLi2]=get_interpolation_on_grid(method,COAST.xq,COAST.yq,COAST.x,COAST.y,SWLi,[]);
        [~,closest]=min(hypot((x_struc(1)+x_struc(2))/2-COAST.xq,(y_struc(1)+y_struc(2))/2-COAST.yq));
        SWLstruc = SWLi2(closest); 

        %% Calculate relevant variables
        kh=wave_GUO2002(Tp,BW_depth);
        Lwave=2*pi*BW_depth./kh;
        Lwave_vGent=(9.81/(2*pi))*((Tp/1.2)^2);
        Rc=BW_cheight-SWLstruc;  % Rc is positive when the crest is higher than the water level, meaning lower Kt for higher Rc. Typically it would be between -2m below water and +1m above water level
        CHIop=BW_slope./((Hs./Lwave).^0.5);
        
        %% Calculate wave transmission
        if ~isempty(findstr(lower(BW_transm_form),'angr'))   % d'Angremond et al. (1996) formulation, rough / permeable structure 
            Kt=real(-0.4.*Rc/Hs+0.64.*((BW_width/Hs).^-0.31).*(1 - exp(-0.5.*CHIop)));
         elseif ~isempty(findstr(lower(BW_transm_form),'gent')) % van Gent et al. (2023) formulation, rough / permeable structure
            c1=0.43; 
            c2=3.1;
            c3=0.75;
            c4=-0.25;
            c5=0.5;
            Kt=c1*tanh(-(Rc/Hs+c2*((BW_width/Lwave_vGent)^c3)+c4))+c5;     
        elseif ~isempty(findstr(lower(BW_transm_form),'seabrhall')) % Seabrook & Hall (1998) formulation, rough / permeable structure 
            Kt=1-exp(-0.65*(-1*Rc/Hs)-1.09*(Hs/BW_width))+0.047*((BW_width*-1*Rc)/(Lwave*BW_D50))-0.067*(-1*Rc*Hs)/(BW_width*BW_D50); 
        end
        Kt=min(Kt,1);   % Limitations are that Kt=0.8 and Kt=0.075 in d'Angremond et al. formula. 
        Kt=max(Kt,0);        
        Kt_total=(Kt_total*(tt-1)+Kt)/tt;   % compute average Kt over various waterlevels 
    end
end

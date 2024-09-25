function [Kt_total]=wave_transmission(xhard,yhard,xq,yq,x,y,SWL,Hs,Tp,transmbwdepth,transmcrestheight,transmcrestwidth,transmslope,transmform,transmd50)
% function [Kt_total]=wave_transmission(xhard,yhard,xq,yq,x,y,SWL,Hs,Tp,transmbwdepth,transmcrestheight,transmcrestwidth,transmslope,transmform,transmd50)
% 
% Wave transmission calculation.
% 
% INPUT: 
%     xhard             : x-coordinates of the offshore breakwater [m]
%     yhard             : y-coordinates of the offshore breakwater [m]
%     xq                : x-coordinates of the transport points [m]
%     yq                : y-coordinates of the transport points [m]
%     x                 : x-coordinates of the coastline points [m]
%     y                 : y-coordinates of the coastline points [m]
%     Hs                : significant wave height [m]
%     Tp                : peak wave period [s]
%     Dir               : wave direction [deg]
%     SWL               : water level [m]
%     transmbwdepth     : depth at breakwater location [m]
%     transmcrestheight : breakwater crest height [m]
%     transmcrestwidth  : breakwater crest width [m]
%     transmslope       : breakwater slope [-]
%     transmform        : used formulation for calculating transmission ('angr', 'gent' or 'seabrhall')
%     transmd50         : D50 of breakwater armour material [m]
% 
% OUTPUT:
%     Kt_total          : transmitted wave [0-1]
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
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses>
%   --------------------------------------------------------------------

    %% Initialize fluxes
    [ x_struc,y_struc,nhard,~,~ ] = get_one_polygon(xhard,yhard,1);
    Kt = zeros(length(Hs),length(nhard));
    Kt_total = zeros(length(Hs),length(nhard));
    
    % Use sub-timesteps when multiple measurements of SWL are taking place in a single timestep 
    % (i.e. at more than just one timestep, but also inbetween moments)
    dtfractions=max([size(Hs,1),size(Tp,1),size(SWL,1)]);
  
    % Loop over the sub-conditions of the time step (i.e. at fractions of coastline timestep)
    for tt=1:dtfractions   
        
        %% Prepare SWL format
        SWLi=SWL(min(tt,size(SWL,1)),:);   
        % Interpolate SWL conditions along coast for wave locations
        method='weighted_distance';
        [SWLi2]=get_interpolation_on_grid(method,xq,yq,x,y,SWLi,[]);
        [~,closest]=min(hypot((x_struc(1)+x_struc(2))/2-xq,(y_struc(1)+y_struc(2))/2-yq));
        SWLstruc = SWLi2(closest); 

        %% Calculate relevant variables
        h=transmbwdepth;
        k=get_disper(h,Tp);
        Lwave=2*pi./k;
        Lwave_vGent=(9.81.*(Tp/1.2).^2)/(2*pi);
        Rc=transmcrestheight-SWLstruc;  % Rc is positive when the crest is higher than the water level, meaning lower Kt for higher Rc. Typically it would be between -2m below water and +1m above water level
        CHIop=transmslope./((Hs./Lwave).^0.5);
        
        %% Calculate wave transmission
        if ~isempty(findstr(lower(transmform),'angr'))   % d'Angremond et al. (1996) formulation, rough / permeable structure 
            Kt=real(-0.4.*Rc/Hs+0.64.*((transmcrestwidth/Hs).^-0.31).*(1 - exp(-0.5.*CHIop)));
         elseif ~isempty(findstr(lower(transmform),'gent')) % van Gent et al. (2023) formulation, rough / permeable structure
            c1=0.43; 
            c2=3.1;
            c3=0.75;
            c4=-0.25;
            c5=0.5;
            Kt=c1*tanh(-(Rc/Hs+c2*((transmcrestwidth/Lwave_vGent)^c3)+c4))+c5;     
        elseif ~isempty(findstr(lower(transmform),'seabrhall')) % Seabrook & Hall (1998) formulation, rough / permeable structure 
            Kt=1-exp(-0.65*(-1*Rc/Hs)-1.09*(Hs/transmcrestwidth))+0.047*((transmcrestwidth*-1*Rc)/(Lwave*transmd50))-0.067*(-1*Rc*Hs)/(transmcrestwidth*transmd50); 
        end
        Kt=min(Kt,1);   % Limitations are that Kt=0.8 and Kt=0.075 in d'Angremond et al. formula. 
        Kt=max(Kt,0);        
        Kt_total=(Kt_total*(tt-1)+Kt)/tt;   % compute average Kt over various waterlevels 
    end
end

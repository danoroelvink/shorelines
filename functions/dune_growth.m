function [qw,uc,us,mwe,qwe,Bdry,FL]=dune_growth(uwindi,duneaw,rhos,rhoa,por,g,d50,d50r,z,k,kw,R,SWLi,dfelev,wberm,phic,phiwndi,segmaw)
% function [qw,uc,us,mwe,qwe,Bdry,FL]=dune_growth(uwindi,duneaw,rhos,rhoa,por,g,d50,d50r,z,k,kw,R,SWLi,dfelev,wberm,phic,phiwndi,segmaw)
% 
% The rate of dune growth due to aeolian transport is computed in this routine.
% 
% INPUT:
%    uwindi         : wind velocity [m/s]
%    duneaw         : Coefficient (Bagnold, 1937)
%    rhos           : Density of the sediment [kg/m3]
%    rhoa           : Density of air [kg/m3]
%    por            : porosity
%    g              : Gravitational acceleration [m/s2]
%    d50            : Median grain diameter [m]
%    d50r           : Median reference grain size [m]
%    z              : wind measurement vertical height [m]
%    k              : Von Karman's coefficient
%    kw             : Empirical coefficient (Sherman et al. 2013)
%    R              : runup level for dune erosion [m]
%    SWLi           : still water level for dune erosion [m]
%    dfelev         : Dune foot elevation [m]
%    wberm          : Berm width (distance MSL to dune foot) [m]
%    phic           : shore-normal orientation at coastline points [°N]
%    phiwndi        : wind direction [°N]
%    segmaw         : Empirical factor used for scaling impact of the fetch length
% 
% OUTPUT:
%    wberm           : berm width (distance MSL to dune foot) [m]
%    qw              : wind transport from the beach to the dune [m3/m/yr]
%    uc              : critical and shear velocity
%    us              : wind shear velocity (m/s)
%    mwe             : potential aeolian transport rate (Kg/m/s)
%    qwe             : equilibrium transport rate [m3/m/yr] 
%    Bdry            : dry beach width [m]
%    FL              : fetch length [m]
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
   
    n=length(phic);
    
    % Critical and shear velocity
    uc=duneaw*sqrt((rhos-rhoa)*g*d50/rhoa);                                     % wind critical shear velocity (m/s)
    uc=uc.*ones(1,n);
    z0=d50/30;                                                                  % Larson 2016
    %z0=0.081*log(DUNE.d50/0.18);                                               % Zingg (1953) Hallin et al. (2019)
    us=k.*uwindi./log(z/z0);                                                    % wind shear velocity (m/s)
    
    % Potential aeolian transport rate (Kg/m/s)
    %mov=us>uc;
    mwe=kw*sqrt(d50r/d50)*rhoa*us.^2.*max(us-uc,0)/g;               
    %mwe(~mov)=0;
    
    % Equilibrium transport rate [m3/m/yr] 
    spyr     = 3600*24*365;
    qwe=spyr*mwe/rhos/(1-por);                                                  % adapted by Anouk                                                 
    
    % Estimation of the dry beach width
    Bdry=max((1-(R+SWLi)./dfelev),0).*wberm;      
    
    % Fetch length
    cosphiwndloc=max(cosd(phic-phiwndi),0);
    FL = max(Bdry./cosphiwndloc,0);
    
    % Aeolian transport rate
    qw0=qwe.*(1-(exp(-segmaw.*FL))); 
    qw=qw0.*cosphiwndloc; 
end
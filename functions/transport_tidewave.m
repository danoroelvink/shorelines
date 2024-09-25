function [Slongtot,Slong_mean,Slong,v,vt,vw,Hrms,h]=transport_tidewave(eta,detads,phi,surfslope,Ttide,nT,ktide,cf,hmin,Hrms0,Tp,theta0,alpha,gamma,ks,D50,D90,hclosure,x,zb,Acal,n)
% function [Slongtot,Slong_mean,Slong,v,vt,vw,Hrms,h]=transport_tidewave(eta,detads,phi,surfslope,Ttide,nT,ktide,cf,hmin,Hrms0,Tp,theta0,alpha,gamma,ks,D50,D90,hclosure,x,zb,Acal,n)
% 
% TIDE_1D computes the tide-wave driven transports
% 
% INPUT: 
%     eta         : M2 and M4 water level amplitudes (m)
%     detads      : longshore gradients of M2 and M4 tidal amplitudes
%     phi         : phase (deg) of M2 and M4 water level components
%     surfslope   : longshore mean surface slope driving residual current
%     Ttide       : duration of a spring-neap cycle (745*60)
%     nT          : number of points in a single tide (24)
%     ktide       : longshore wave number M2 and M4 (rad/m)
%     cf          : roughness factor [-]
%     hmin        : minimum depth [m]
%     Hrms0       : Square root wave height at nearshore location (or diffracted wave height) [Nx1]
%     Tp          : Wave period [s]
%     theta0      : Relative angle of waves with respect to the coastline at the depth-of-closure [°]
%     alpha       : calibration factor for point of braking [-]
%     gamma       : depth-induced breaking coefficient [-]
%     ks          : roughness height [m] 
%     d50         : Median grain size [um]
%     d90         : 90-th percentile grain size [um]
%     hclosure    : depth-of-closure [m]
%     x           : x-coordinates of cross-shore profiles [m]
%     zb          : z-coordinates of cross-shore profiles [m]
% 
% OUTPUT:
%     Slongtot    : total longshore transport [m3/yr]
%     Slong_mean  : longshore transport over the cross-shore profile averaged over the tide [nx] [m3/m/yr]
%     Slong       : longshore transport over the cross-shore profile for all tide timesteps [nx nt] [m3/m/yr]
%     v           : alongshore velocity [nx nT] [m/s]
%     vt          : alongshore velocity due to the tide [nx nT] [m/s]
%     vw          : longshore wave-driven current velocity, 2D matrix [nx nT] [m/s]
%     Hrms        : Hrms wave height, 2D matrix [nx nT] (m)
%     h           : depth over the cross-shore profile [nx nT] [m]
% 
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2024 IHE Delft & Deltares
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

    %% Compute tidal velocity, intra-tide and residual
    [h,vt] = tide_1d_ana_anycomp(eta,detads,phi,0,surfslope,Ttide,nT,ktide,cf,hmin,x,zb);
    
    %% Compute wave decay throughout the tide and wave-driven current 
    [Hrms,Dw,Urms,k,C,Cg,theta,Fy,vw] = wave_cur_1D(Hrms0,Tp,theta0,alpha,gamma,cf,x,h,n);
    
    %% Total current (x,t)
    v=vt+vw;
    
    %% Longshore transport for each x and t
    [Slong]=transport_soulsbyvanrijn(h,Tp,k,Hrms,v,ks,hmin,D50,D90,Acal);
    
    %% Average longshore transport profile
    Slong(h==hmin)=0;
    Slong(h>hclosure)=0;
    Slong_mean=mean(Slong,2);
    hmean=mean(h,2);
    
    %% Total longshore transport
    Slongtot=sum(Slong_mean)*(x(2)-x(1));
    
end


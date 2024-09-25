function [Slong]=transport_soulsbyvanrijn(h,T,k,H,v,ks,hmin,D50,D90,Acal)
% function [Slong]=transport_soulsbyvanrijn(h,T,k,H,v,ks,hmin,D50,D90,Acal)
% 
% Compute sand transport according to Soulsby - van Rijn
% 
% INPUT: 
%     h           : depth over the cross-shore profile [nx nT] [m]
%     T           : Wave period [s]
%     k           : wave number, 2D matrix [nx nT] (1/m)
%     H           : Hrms wave height, 2D matrix [nx nT] (m)
%     v           : alongshore velocity [nx nT] [m/s]
%     ks          : roughness height [m] 
%     hmin        : minimum depth [m]
%     d50         : Median grain size [um]
%     d90         : 90-th percentile grain size [um]
% 
% OUTPUT:
%     Slong       : longshore transport over the cross-shore profile for all tide timesteps [nx nT] [m3/m/yr]
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

    g=9.81;delta=1.65;rho=1025;nu=1.e-6;kappa=0.39;
    z0=ks/30;
    %z0=0.006;
    Dstar=(g*delta/nu^2)^(1/3)*D50;
    dry=h<hmin;
    h=max(h,hmin);
    Urms=1/sqrt(2)*pi*H./T./sinh(k.*h);
    if D50<=0.5e-3
        Ucr=0.19*D50^0.1*log10(4*h/D90);
    else
        Ucr=8.5*D50^0.6*log10(4*h/D90);
    end
    cf=(kappa./(log(h/z0)-1)).^2;
    umod=sqrt(v.^2+0.018./cf.*Urms.^2);
    ksi=(umod-Ucr).^2.4;
    ksi(umod<Ucr)=0;
    Asb=0.005*h.*(D50./h/delta/g/D50).^1.2;
    Ass=0.012*D50*Dstar^(-0.6)/(delta*g*D50)^1.2;
    Sby=Acal*Asb.*v.*ksi;
    Ssy=Acal*Ass.*v.*ksi;
    Stoty=Sby+Ssy;
    %% Set values to 0 at dry points
    Urms(dry)=0;
    Ucr(dry)=0;
    Cd(dry)=0;
    ksi(dry)=0;
    Asb(dry)=0;
    Ass(dry)=0;
    Sbx(dry)=0;
    Sby(dry)=0;
    Ssx(dry)=0;
    Ssy(dry)=0;
    Stotx(dry)=0;
    Stoty(dry)=0;
    Slong=Stoty*3600*24*365;
end

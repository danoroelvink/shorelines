function [TRANSP,COAST] = transport_mud(COAST,TRANSP,WAVE,WIND,MUD)
% function [TRANSP,COAST] = transport_mud(COAST,TRANSP,WAVE,WIND,MUD)
% 
% This function computes the transport rate of mud along the coast,
% on the basis of alongshore discharges of water due to tide, waves and wind.
% Concentrations are computed based on alongshore supply and erosion
% of the cohesive bed. The coastline change is computed, as well as 
% the rate of change of the mud flats and mangroves. 
% 
% INPUT:
%    COAST
%       .PHIc           : angle of coastline [°N]
%       .ds             : grid cell size [m]
%       .h0             : active height [m]
%       .Bf             : mudflat width [m]
%       .Bm             : mangrove width [m]
%       .Bfm            : colonizing mangrove width [m]
%    TRANSP
%       .rhow           : density of water [kg/m3]
%       .rhos           : density of sediment [kg/m3]
%       .cf             : roughness coefficient [m]
%       .g              : acceleration of gravity [m/s2]
%    WAVE
%       .HStdp          : significant wave height at depth-of-closure [m]
%       .dPHItdp        : relative wave direction w.r.t. coastline at depth-of-closure [°]
%       .dnearshore     : nearshore depth (depth-of-closure) [m]
%    WIND
%       .phiwnd         : incoming direction of the wind [°N]
%       .uz             : velocity of the wind [m/s]
%       .rhoa           : density of air [kg/m3]
%       .Cd             : drag coefficient of wind force on water [-]
%    MUD
%       .used           : switch for mud transport option (0/1)
%       .taucr          : critical shear stress for erosion [N/m2]
%       .B              : width of the mud transport zone [m]
%       .M              : erosion rate [kg/m2/s]
%       .w              : fall velocity of the muddy sediment [m/s] 
%       .Bfcrit         : critical mudflat width [m]
%       .MHW            : mhw level [m]
%       .MSL            : msl level [m]
%       .dm             : depth at mud flat interface [m]
%       .Tfm            : TFM parameter
%       .rate_density   : riverdischarge rate for each individual measure in [m3/yr]
%
% OUTPUT:
%     TRANSP    
%       .C              : concentration of mud in the water [kg/m3]
%       .Q              : flow of water along the coast [m3/s]
%       .QS             : transport rate [m3/yr]
%     COAST    
%       .dndt_mud       : coastline change rate due to mud transport gradients [m/yr]
%       .dBfdt          : change in mudflat width [m/yr]
%       .dBfmdt         : change in mangove width [m/yr]
%       .dBmdt          : change in colonizing mangrove width [m/yr]
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

    if MUD.used==1
        phic     = COAST.PHIc;
        ds       = COAST.ds;
        d        = COAST.h0;
        nq       = length(phic);
        n        = nq-1;
        phiwnd   = WIND.phiwnd;
        uwind    = WIND.uz;
        rhoa     = WIND.rhoa;
        Cd       = WIND.Cd;
        rhow     = TRANSP.rhow;
        rhos     = TRANSP.rhos;
        cf       = TRANSP.cf;
        g        = TRANSP.g;
        taucr    = MUD.taucr;
        B        = MUD.B;
        M        = MUD.M;
        w        = MUD.w;
        Bfcrit   = MUD.Bfcrit;
        MHW      = MUD.MHW;
        MSL      = MUD.MSL;
        Tfm      = MUD.Tfm;
        dm       = MUD.dm;
        qriv     = MUD.rate_density;
        spyr     = 3600*24*365;

        %% Mangrove properties
        Bf       = COAST.Bf;     % Mudflat width
        Bm       = COAST.Bm;     % Mangrove width
        Bfm      = COAST.Bfm;    % Colonizing mangrove width
      
        %% Computing tidal prism and exchange
        P        = (0.5*Bm+Bfm)*(MHW-MSL);         % Tidal prism in mangrove
        qmoverC  = 700/spyr*P.*max(min((Bf-Bfcrit)/Bfcrit,1),0);
        % transport into mangrove / C

        %% Waves and current
        WAVE.HStdp(abs(WAVE.dPHItdp)>90) = 0.1;
        urms     = sqrt(g/8/WAVE.dnearshore)*WAVE.HStdp;
        A        = 0.5*B.*WAVE.dnearshore;
        Fws      = rhow*g/32*WAVE.HStdp.^2.*sind(2*(WAVE.dPHItdp))/B;
        arg      = (rhoa*Cd*uwind.^2.*sind(phic-phiwnd)+Fws)/rhow/cf./urms;
        aa       = 1./urms.^2;bb=1.16^2;cc=-arg.^2;
        v        = sqrt((-bb+sqrt(bb.^2-4*aa.*cc))/2./aa).*sign(arg);
        Q        = A.*v;
        for i=1:n
           im1=max(i-1,1);
           ip1=min(i+1,nq);
           if Q(i)>0 && Q(ip1)>0
              a(i)=-Q(i)/ds(i);
              b(i)=Q(ip1)/ds(i)+w*B+qmoverC(i);
              c(i)=0;
           elseif Q(i)<0 && Q(ip1)<0
              a(i)=0;
              b(i)=-Q(i)/ds(i)+w*B+qmoverC(i);
              c(i)=Q(ip1)/ds(i);
           elseif Q(i)>0 && Q(ip1)<0
              a(i)=-Q(i)/ds(i);
              b(i)=w*B+qmoverC(i);
              c(i)=Q(ip1)/ds(i);
           else
              a(i)=0;
              b(i)=(Q(ip1)-Q(i))/ds(i)+w*B+qmoverC(i);
              c(i)=0;
           end
           vm=abs(.5*(v(i)+v(ip1)));
           %tau(i)=cf*rhow*vm*max(vm,urms(i));
           tau(i)=cf*rhow*(max(vm,urms(i))^2);
           r(i)=M*B*max((tau(i)-taucr),0)/taucr+qriv(i)*rhos/spyr;
        end
        a(1)=0;b(1)=-1;c(1)=1;r(1)=0;
        a(n)=-1;b(n)=1;c(n)=0;r(n)=0;
        C=tridiag(a,b,c,r,n);
        if min(C)<0
           disp('negative concentrations')
        end
        C=max(C,0);
        qm=qmoverC.*C;
        TRANSP.Q  = Q; % discharges along the coast [m3/s]
        TRANSP.C  = C; % concentrations [kg/m3]
        TRANSP.QS = zeros(size(Q));
        for i=2:length(TRANSP.Q)-1
           if TRANSP.Q(i)>0
              TRANSP.QS(i)=(TRANSP.C(i-1)*TRANSP.Q(i))/rhos*spyr; % [m3/yr]
           else
              TRANSP.QS(i)=(TRANSP.C(i)*TRANSP.Q(i))/rhos*spyr;   % [m3/yr]
           end
        end
        TRANSP.QS(1)=(TRANSP.Q(1)*TRANSP.C(1))/rhos*spyr;         % [m3/yr]
        TRANSP.QS(end)=(TRANSP.Q(end)*TRANSP.C(end))/rhos*spyr;   % [m3/yr]
        dQsds=diff(TRANSP.QS)./ds;                                % [m3/m/yr]
        qm=qm/rhos*spyr;                                          % [m3/m/yr]
        for i=1:n
           if Bf(i)>0 || -dQsds(i)-qm(i)+qriv(i)>0
              dBmdt(i)=2*qm(i)/dm;
              dndt(i)=(-dQsds(i)-qm(i)+qriv(i))/d(i);
              dBfdt(i)=dndt(i)-dBmdt(i);
              dBfmdt(i)=(Bf(i)-Bfm(i))/Tfm;
           elseif Bf(i)==0 && -dQsds(i)-qm(i)+qriv(i)<=0
              if Bfm(i)==0
                 dndt(i)=(-dQsds(i)+qriv(i))/(d(i)+0.5*dm);
                 dBmdt(i)=dndt(i);
                 dBfdt(i)=0;
              else
                 dndt(i)=-dQsds(i)/d(i);
                 dBmdt(i)=0;
                 dBfdt(i)=0;
                 dBfmdt(i)=dndt(i);
              end
           end
        end
        COAST.dndt_mud = dndt;
        COAST.dBfdt    = dBfdt;
        COAST.dBfmdt   = dBfmdt;
        COAST.dBmdt    = dBmdt;
    end
end

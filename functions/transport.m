function [TRANSP]=transport(TRANSP,WAVE,TIDE,STRUC,useQSmax)
% function [TRANSP]=transport(TRANSP,WAVE,TIDE,STRUC,useQSmax)
%
%
% INPUT :
%      TRANSP          : Structure with input information of the ShorelineS model, with relevant fields:
%         .dnearshore  : depth at toe of dynamic profile
%         .trform      : transport formulation (either 'CERC', 'KAMP', 'MILH', 'CERC3', 'VR14')
%         .b           : Calbration factor of CERC formula (only CERC)
%         .TP        : Wave period [s]
%         .d50         : Median grain size [um]
%         .tanbeta     : Slope [1:n] (only Kamphuis and Milhomens)
%         .Pswell      : Relative part of the Pswell (in percentage, specify a value between 0 and 100)
%         .rhos        : Density of sediment [kg/m3]
%         .rhow        : Density of water [kg/m3]
%         .gamma       : Breaking coefficient [-]
%         .porosity    : Porosity [-]
%         .g           : Acceleration of gravity [m/s2]
%      WAVE
%         .dPHItdp     : Relative angle of waves with respect to the coastline (corrected to 0 in shadow zones) at nearshore location (or offshore for CERC1 and CERC2)
%         .dPHIbr      : Relative angle of waves with respect to the coastline (corrected to 0 in shadow zones) at the point of breaking
%         .HStdp       : Wave height at nearshore location (or diffracted wave height) [Nx1]
%         .HSbr        : Breaking wave height (or diffracted wave height) [Nx1]
%         .c1          : inclination dQS/dTHETA of s-phi curve (in case of rays)
%         .c2          : curvature coefficient of s-phi curve (in case of rays)
%         .QSoffset    : transport offset of the s-phi curve (in case of rays)
%         .diffraction : switch for wave diffraction
%      TIDE
%         .eta         : M2 and M4 water level amplitudes (m)
%         .detads      : longshore gradients of M2 and M4 tidal amplitudes
%         .phi         : phase (deg) of M2 and M4 water level components
%         .k           : longshore wave number M2 and M4 (rad/m)
%         .ss          : longshore mean surface slope driving residual
%                        current
%      STRUC
%         .xp          : x-coordinate of diffraction point
%         .yp          : y-coordinate of diffraction point
%      useQSmax        : switch for using QSmax (either use '1' or 'QSmax' or 'yes' as input to switch output to field '.QSmax' on) 
% 
% OUTPUT :
%      TRANSP
%         .QS          : Transport rates in grid cells [1xN] (in [m3/yr] including pores)
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
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses
%   --------------------------------------------------------------------
    %fprintf('  Transport computation \n');

    fieldQS='QS';
    if nargin>4
        if isnumeric(useQSmax)
            if useQSmax~=0
                fieldQS='QSmax';
            end
        else
            if ~strcmpi(useQSmax,'no')
                fieldQS='QSmax';
            end
        end
    end
    %% Initialize QS
    QS=zeros(size(WAVE.HStdp));
    
    %% Compute gradient of HS
    nq=length(WAVE.HStdp);
    dHS=zeros(1,nq);
    if WAVE.diffraction==1 && strcmpi(fieldQS,'QS')
        HS=WAVE.HStdp;
        if ~strcmpi(TRANSP.trform,'CERC') && ~strcmpi(TRANSP.trform,'CERC2') && ~strcmpi(TRANSP.trform,'KAMP') && ~strcmpi(TRANSP.trform,'MILH')
            HS=WAVE.HSbr;
        end
        dHS=zeros(size(HS));
        xq = STRUC.xq;
        yq = STRUC.yq;      
        if nq>1
            for i=1:nq
                im1=max(i-1,1);
                ip1=min(i+1,nq);
                dHS(i)=(HS(ip1)-HS(im1)) / hypot(yq(ip1)-yq(im1),xq(ip1)-xq(im1)) ;
            end
        end
    end 
    
    %% Transport : CERC with the offshore wave height and direction
    if strcmpi(TRANSP.trform,'CERC')
        k=0.2;                                                       % using CERC (1984) value of k(SPM,Hs)=0.39 is suggested, but this is typically quite high
        QS = TRANSP.qscal .* TRANSP.b .* WAVE.HStdp.^2.5 .* (sind(2*WAVE.dPHItdp)-2*cosd(WAVE.dPHItdp).* dHS); % use HS and dPHI, and dHS as second-order component dHS   
        QS(abs(WAVE.dPHItdp)>90)=0;     % Set transport to 0 when angle exceeds 180 degrees
    
    %% Transport : CERC with the offshore wave height and direction, including an implicit refraction from offshore to nearshore within the transport formulation
    elseif strcmpi(TRANSP.trform,'CERC2')
        k=0.39;                                                      % using CERC (1984) value of k(SPM,Hs)=0.39
        b1 = k .* (TRANSP.rhow .* TRANSP.g.^0.5 ./ (16 .* sqrt(TRANSP.gamma).* (TRANSP.rhos-TRANSP.rhow) .* (1-TRANSP.porosity)));    % = theoretical 'b1' factor from Shore Protection Manual
        b2 = b1 .* ((TRANSP.gamma.*TRANSP.g).^0.5 ./(2*pi)).^0.2;               % b_theory = 0.0946 * 365*24*60*60 = 2.9833E+6
        QS = TRANSP.qscal .* 365*24*60*60*b2.*WAVE.HStdp.^(12/5).*WAVE.TP.^(1/5).* (cosd(WAVE.dPHItdp).^(6/5).*sind(WAVE.dPHItdp)-2*cosd(WAVE.dPHItdp).* dHS); % use HS and dPHI, and dHS as second-order component dHS
        QS(abs(WAVE.dPHItdp)>90)=0;     % Set transport to 0 when angle exceeds 180 degrees
    
    %% Transport: CERC with the nearshore breaking wave height and direction (e.g. Vitousek, Barnard, 2015)
    elseif strcmpi(TRANSP.trform,'CERC3')
        k=0.35;                                                                                  % using a default value of k=0.35
        b3 = 1/16 * k  .* TRANSP.rhow ./ ((TRANSP.rhos-TRANSP.rhow) .* (1-TRANSP.porosity));     % b3 by default is 0.023 for CERC3 using k=0.35, rhos=2650, rhow=1025 and porosity=0.4
        QS = TRANSP.qscal .* 365*24*60*60 .* b3 .* TRANSP.g.^0.5 .* (TRANSP.gamma).^-0.52 .* WAVE.HSbr.^2.5 .* (sind(2.*WAVE.dPHIbr)-2*cosd(WAVE.dPHIbr).* dHS); % use HS and dPHI, and dHS as second-order component dHS
        QS(abs(WAVE.dPHIbr)>90)=0;      % Set transport to 0 when angle exceeds 180 degrees
    
    %% Transport : Kamphuis
    elseif strcmpi(TRANSP.trform,'KAMP')
        % use HS and dPHI, and dHS as second-order component
        %QSkampmass=TRANSP.qscal .* 2.33 .* WAVE.TP.^1.5 .* TRANSP.tanbeta.^0.75 .* TRANSP.d50.^-0.25 .* WAVE.HStdp.^2 .* ( abs(sind(2*WAVE.dPHIbr)).^0.6.*sign(WAVE.dPHIbr));
        %QSkampmass=TRANSP.qscal .* 2.33 .* WAVE.TP.^1.5 .* TRANSP.tanbeta.^0.75 .* TRANSP.d50.^-0.25 .* WAVE.HStdp.^2 .* ( abs(sind(2*WAVE.dPHIbr)).^0.6.*sign(WAVE.dPHIbr) - 2.*cosd(WAVE.dPHIbr).*dHS );
        QSkampmass=TRANSP.qscal .* 2.33 .* WAVE.TP.^1.5 .* TRANSP.tanbeta.^0.75 .* TRANSP.d50.^-0.25 .* WAVE.HSbr.^2 .* ( abs(sind(2*WAVE.dPHIbr)).^0.6.*sign(WAVE.dPHIbr) - (2/TRANSP.tanbeta).*cosd(WAVE.dPHIbr).*dHS );
        QS = TRANSP.qscal .* 365*24*60*60*(QSkampmass /(TRANSP.rhos-TRANSP.rhow)) /(1.0-TRANSP.porosity);
        QS(abs(WAVE.dPHIbr)>90)=0;      % Set transport to 0 when angle exceeds 180 degrees
    
    %% Transport : Mil-Homens (2013)
    elseif strcmpi(TRANSP.trform,'MILH')
        QSmilhmass=TRANSP.qscal .* 0.15.* WAVE.TP.^0.89 .* TRANSP.tanbeta.^0.86 .* TRANSP.d50.^-0.69 .* WAVE.HStdp.^2.75 .*( abs(sind(2*WAVE.dPHIbr)).^0.5.*sign(WAVE.dPHIbr) - 2.*cosd(WAVE.dPHIbr).*dHS ); % <- this is the immersed mass under water
        QS = TRANSP.qscal .* 365*24*60*60*(QSkampmass /(TRANSP.rhos-TRANSP.rhow)) /(1.0-TRANSP.porosity);
        QS(abs(WAVE.dPHIbr)>90)=0;     % Set transport to 0 when angle exceeds 180 degrees
    
    %% Transport : Van Rijn (2014)        
    elseif strcmpi(TRANSP.trform,'VR14')
        vtide=0;
        kswell=0.015.*TRANSP.Pswell+(1-0.01.*TRANSP.Pswell);
        vwave=0.3.*real(sind(2.*WAVE.dPHIbr)).*(TRANSP.g.*WAVE.HSbr).^0.5; % added secondary component due to wave height difference
        %vwave=0.3.*real(sind(2.*WAVE.dPHIbr)-2.*cosd(WAVE.dPHIbr).*dHS).*(TRANSP.g.*WAVE.HSbr).^0.5; % added secondary component due to wave height difference
        %vwave=0.3.*real(sind(2.*WAVE.dPHIbr)-2/TRANSP.tanbeta.*cosd(WAVE.dPHIbr).*dHS).*(TRANSP.g.*WAVE.HSbr).^0.5; % added secondary component due to wave height difference
        vtotal=vwave+vtide;
        QSvr14mass=TRANSP.qscal .* 0.0006 .* kswell .* TRANSP.rhos .* TRANSP.tanbeta.^0.4 .*TRANSP.d50.^-0.6 .* WAVE.HSbr.^2.6 .* abs(vtotal).*sign(WAVE.dPHIbr);
        QS = 365*24*60*60*(QSvr14mass ./(TRANSP.rhos-TRANSP.rhow)) ./(1.0-TRANSP.porosity); 
        QS(abs(WAVE.dPHIbr)>90)=0;     % Set transport to 0 when angle exceeds 180 degrees
    
    elseif strcmpi(TRANSP.trform,'RAY')
        QS = -WAVE.c1.*WAVE.dPHItdp .* exp(-(WAVE.c2.*WAVE.dPHItdp).^2) + WAVE.QSoffset;
        QS(abs(WAVE.dPHItdp)>90)=0;     % Set transport to 0 when angle exceeds 180 degrees
        dQSdPHI = WAVE.c1.*exp(-(WAVE.c2.*WAVE.dPHItdp).^2) .* ( -2.*(WAVE.c2.*WAVE.dPHItdp).^2 + 1. ) *180./pi; 

    elseif strcmpi(TRANSP.trform,'WAVETIDEPROF')
        alpha=1.0;    % Default dissipation parameter Baldock
        gamma=0.78;   % Default breaking parameter Baldock
        for i=1:length(WAVE.HStdp)
            Hrms0=WAVE.HStdp(i)/sqrt(2);
            Tp=WAVE.TP(i);
            theta0=WAVE.dPHItdp(i);
            [QS(i),Slong_mean,Slong,v,vt,vw,Hrms,h] = tide_wave_transport ...
            (TIDE.eta(i,:),TIDE.detads(i,:),TIDE.phi(i,:),TIDE.ss(i), ...
             TIDE.Ttide,TIDE.nT,TIDE.k(i,:),TIDE.Cf,TIDE.hmin, ...
             Hrms0,Tp,theta0,alpha,gamma, ...
             TRANSP.ks,TRANSP.d50,TRANSP.d90,TIDE.hclosure, ...
             TIDE.x,TIDE.zb);
        end

    end

    %% Set transport to 0 when angle exceeds 180 degrees
    % QS(abs(dPHI)>90)=0; 

    TRANSP.(fieldQS)=QS;
    if strcmpi(fieldQS,'QS')
        TRANSP.debug.QS0=TRANSP.QS;
    end
end

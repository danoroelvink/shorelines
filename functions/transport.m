function [TRANSP]=transport(TRANSP,WAVE,TIDE,STRUC,useQSmax)
% function [TRANSP]=transport(TRANSP,WAVE,TIDE,STRUC,useQSmax)
%
% This function computes the alongshore wave or tide driven transport rates. 
% Multiple transport formulae can be used. 
%
% INPUT: 
%      TRANSP          : Structure with input information of the ShorelineS model, with relevant fields:
%         .dnearshore  : depth at toe of dynamic profile [m]
%         .trform      : transport formulation (either 'CERC', 'KAMP', 'MILH', 'CERC3', 'VR14')
%         .qscal       : calibration coefficient for the transport computations [-]
%         .b           : transport factor of CERC formula (only CERC)
%         .d50         : median grain size [um]
%         .d90         : 90-th percentile grain size [um]
%         .ks          : roughness height [m] 
%         .tanbeta     : slope [1:n] (only Kamphuis and Milhomens)
%         .pswell      : relative part of the pswell (in percentage, specify a value between 0 and 100)
%         .rhos        : density of sediment [kg/m3]
%         .rhow        : density of water [kg/m3]
%         .gamma       : depth-induced breaking coefficient [-]
%         .alpha       : calibration factor for point of braking [-]
%         .porosity    : porosity [-]
%         .g           : acceleration of gravity [m/s2]
%      WAVE
%         .HStdp       : wave height at nearshore location (or diffracted wave height) [Nx1]
%         .HSbr        : breaking wave height (or diffracted wave height) [Nx1]
%         .TP          : wave period [s]
%         .dPHItdp     : relative angle of waves with respect to the coastline (corrected to 0 in shadow zones) at nearshore location (or offshore for CERC1 and CERC2)
%         .dPHIbr      : relative angle of waves with respect to the coastline (corrected to 0 in shadow zones) at the point of breaking
%         .c1          : inclination dQS/dTHETA of s-phi curve (in case of rays)
%         .c2          : curvature coefficient of s-phi curve (in case of rays)
%         .QSoffset    : transport offset of the s-phi curve (in case of rays)
%         .diffraction : switch for wave diffraction
%      TIDE
%         .eta         : M2 and M4 water level amplitudes (m)
%         .detads      : longshore gradients of M2 and M4 tidal amplitudes
%         .phi         : phase (deg) of M2 and M4 water level components
%         .k           : longshore wave number M2 and M4 (rad/m)
%         .ss          : longshore mean surface slope driving residual current
%         .Ttide       : duration of a spring-neap cycle (745*60)
%         .nT          : number of points in a single tide (24)
%         .cf          : roughness factor [-]
%         .hmin        : minimum depth [m]
%         .hclosure    : depth-of-closure [m]
%         .x           : x-coordinates of cross-shore profiles [m]
%         .zb          : z-coordinates of cross-shore profiles [m]
%      STRUC
%         .xp          : x-coordinate of diffraction point
%         .yp          : y-coordinate of diffraction point
%      useQSmax        : switch for using QSmax (either use '1' or 'QSmax' or 'yes' as input to switch output to field '.QSmax' on) 
% 
% OUTPUT: 
%      TRANSP
%         .QS          : transport rates in grid cells [1xN] (in [m3/yr] including pores)
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
    nw=size(WAVE.HStdp,1); % number of wave conditions taken along at this timestep (can be more than 1 in case of simultaneous wave conditions)
    nq=size(WAVE.HStdp,2); % number of coastline points
    QS=zeros(nw,nq);
    
    %% Compute gradient of HS
    dHS=zeros(nw,nq);
    if WAVE.diffraction==1 && strcmpi(fieldQS,'QS')
        HS=WAVE.HStdp;
        if ~strcmpi(TRANSP.trform,'CERC') && ~strcmpi(TRANSP.trform,'CERC2') && ~strcmpi(TRANSP.trform,'KAMP') && ~strcmpi(TRANSP.trform,'MILH')
            HS=WAVE.HSbr;
        end
        dHS=zeros(nw,nq);
        xq = STRUC.xq;
        yq = STRUC.yq;      
        if nq>1 || nw>1
            for i=1:nq
                im1=max(i-1,1);
                ip1=min(i+1,nq);
                for j=1:nw
                    dHS(j,i)=(HS(j,ip1)-HS(j,im1)) / max(hypot(yq(ip1)-yq(im1),xq(ip1)-xq(im1)),1) ;
                end
            end
        end
        dHS(isnan(dHS))=0;
    end 
    
    %% Set wave angle to the critical wave angle for transport in case of high-angle incidence
    if strcmpi(fieldQS,'QS') 
        if TRANSP.suppresshighangle==1 
            id = find(WAVE.dPHItdp>WAVE.dPHIcrit);
            dHS(id)=dHS(id).*0;
            dPHItdp=WAVE.dPHItdp;
            dPHItdp(dPHItdp<-90)=-90;
            dPHItdp(dPHItdp>90)=90;
            WAVE.HStdp(abs(dPHItdp)>90)=0;
            %WAVE.HStdp(id)=WAVE.HStdp(id).*max((1+min((abs(dPHItdp(id))-WAVE.dPHIcrit(id))./(WAVE.dPHIcrit(id)-90),0)),0); % linearly decreas the wave energy from dPHIcrit to 90°
            WAVE.dPHItdp(id)=WAVE.dPHIcrit(id);
            WAVE.dPHIbr(id)=WAVE.dPHIcritbr(id);
        end
    end
    
    %% Transport : CERC with the offshore wave height and direction
    if strcmpi(upper(TRANSP.trform),'CERC') || strcmpi(upper(TRANSP.trform),'CERC1')
        k=0.2;                                                       % using CERC (1984) value of k(SPM,Hs)=0.39 is suggested, but this is typically quite high
        %QS = TRANSP.qscal .* TRANSP.b .* WAVE.HStdp.^2.5 .* (sind(2*WAVE.dPHItdp)-2*cosd(WAVE.dPHItdp).* dHS); % use HS and dPHI, and dHS as second-order component dHS (note that the HStdp here is equal to HSoffshore, with the addition of wave diffraction)   
        QS = TRANSP.qscal .* TRANSP.b .* WAVE.HStdp.^2.5 .* (sind(2*WAVE.dPHItdp)-(2/TRANSP.tanbetasetup)*cosd(WAVE.dPHItdp).* dHS); % use HS and dPHI, and dHS as second-order component dHS (note that the HStdp here is equal to HSoffshore, with the addition of wave diffraction)   
        QS(abs(WAVE.dPHItdp)>90)=0;     % Set transport to 0 when angle exceeds 180 degrees
    
    %% Transport : CERC with the offshore wave height and direction, including an implicit refraction from offshore to nearshore within the transport formulation
    elseif ~isempty(findstr(lower(TRANSP.trform),'cerc2'))
        k=0.39;                                                      % using CERC (1984) value of k(SPM,Hs)=0.39
        b1 = k .* (TRANSP.rhow .* TRANSP.g.^0.5 ./ (16 .* sqrt(TRANSP.gamma).* (TRANSP.rhos-TRANSP.rhow) .* (1-TRANSP.porosity)));    % = theoretical 'b1' factor from Shore Protection Manual
        b2 = b1 .* ((TRANSP.gamma.*TRANSP.g).^0.5 ./(2*pi)).^0.2;               % b_theory = 0.0946 * 365*24*60*60 = 2.9833E+6
        %QS = TRANSP.qscal .* 365*24*60*60*b2.*WAVE.HStdp.^(12/5).*WAVE.TP.^(1/5).* (cosd(WAVE.dPHItdp).^(6/5).*sind(WAVE.dPHItdp)-2*cosd(WAVE.dPHItdp).* dHS); % use HS and dPHI, and dHS as second-order component dHS
        QS = TRANSP.qscal .* 365*24*60*60*b2.*WAVE.HStdp.^(12/5).*WAVE.TP.^(1/5).* (cosd(WAVE.dPHItdp).^(6/5).*sind(WAVE.dPHItdp)-(2/TRANSP.tanbetasetup)*cosd(WAVE.dPHItdp).* dHS); % use HS and dPHI, and dHS as second-order component dHS
        QS(abs(WAVE.dPHItdp)>90)=0;     % Set transport to 0 when angle exceeds 180 degrees
    
    %% Transport: CERC with the nearshore breaking wave height and direction (e.g. Vitousek, Barnard, 2015)
    elseif ~isempty(findstr(lower(TRANSP.trform),'cerc3'))
        k=0.35;                                                      % using a default value of k=0.35
        b3 = 1/16 * k  .* TRANSP.rhow ./ ((TRANSP.rhos-TRANSP.rhow) .* (1-TRANSP.porosity));     % b3 by default is 0.023 for CERC3 using k=0.35, rhos=2650, rhow=1025 and porosity=0.4
        %QS = TRANSP.qscal .* 365*24*60*60 .* b3 .* TRANSP.g.^0.5 .* (TRANSP.gamma).^-0.52 .* WAVE.HSbr.^2.5 .* (sind(2.*WAVE.dPHIbr)-2*cosd(WAVE.dPHIbr).* dHS); % use HS and dPHI, and dHS as second-order component dHS
        QS = TRANSP.qscal .* 365*24*60*60 .* b3 .* TRANSP.g.^0.5 .* (TRANSP.gamma).^-0.52 .* WAVE.HSbr.^2.5 .* (sind(2.*WAVE.dPHIbr)-(2/TRANSP.tanbetasetup)*cosd(WAVE.dPHIbr).* dHS); % use HS and dPHI, and dHS as second-order component dHS
        QScerc3mass = QS .* (1.0-TRANSP.porosity) .* TRANSP.rhos / (365*24*60*60); % kg/s
        QS(abs(WAVE.dPHIbr)>90)=0;      % Set transport to 0 when angle exceeds 180 degrees
    
    %% Transport : Kamphuis
    elseif ~isempty(findstr(lower(TRANSP.trform),'kamp'))
        % use HS and dPHI, and dHS as second-order component
        %QSkampmass=TRANSP.qscal .* 2.33 .* TRANSP.rhos/(TRANSP.rhos-TRANSP.rhow) .* WAVE.TP.^1.5 .* TRANSP.tanbeta.^0.75 .* TRANSP.d50.^-0.25 .* WAVE.HStdp.^2 .* ( abs(sind(2*WAVE.dPHIbr)).^0.6.*sign(WAVE.dPHIbr));
        %QSkampmass=TRANSP.qscal .* 2.33 .* TRANSP.rhos/(TRANSP.rhos-TRANSP.rhow) .* WAVE.TP.^1.5 .* TRANSP.tanbeta.^0.75 .* TRANSP.d50.^-0.25 .* WAVE.HStdp.^2 .* ( abs(sind(2*WAVE.dPHIbr)).^0.6.*sign(WAVE.dPHIbr) - 2.*cosd(WAVE.dPHIbr).*dHS ); % added secondary component due to spatial wave height differences 
        QSkampmass=TRANSP.qscal .* 2.33 .* TRANSP.rhos/(TRANSP.rhos-TRANSP.rhow) .* WAVE.TP.^1.5 .* TRANSP.tanbeta.^0.75 .* TRANSP.d50.^-0.25 .* WAVE.HSbr.^2 .* ( abs(sind(2*WAVE.dPHIbr)).^0.6.*sign(WAVE.dPHIbr) - (2/TRANSP.tanbetasetup).*cosd(WAVE.dPHIbr).*dHS ); % added secondary component due to spatial wave height differences + the slope effect on this current
        QS = 365*24*60*60*(QSkampmass / TRANSP.rhos) /(1.0-TRANSP.porosity);
        QS(abs(WAVE.dPHIbr)>90)=0;      % Set transport to 0 when angle exceeds 180 degrees
    
    %% Transport : Mil-Homens (2013)
    elseif ~isempty(findstr(lower(TRANSP.trform),'milh'))
        %QSmilhmass=TRANSP.qscal .* 0.15.*  TRANSP.rhos/(TRANSP.rhos-TRANSP.rhow) .* WAVE.TP.^0.89 .* TRANSP.tanbeta.^0.86 .* TRANSP.d50.^-0.69 .* WAVE.HSbr.^2.75 .* abs(sind(2*WAVE.dPHIbr)).^0.5;
        %QSmilhmass=TRANSP.qscal .* 0.15.*  TRANSP.rhos/(TRANSP.rhos-TRANSP.rhow) .* WAVE.TP.^0.89 .* TRANSP.tanbeta.^0.86 .* TRANSP.d50.^-0.69 .* WAVE.HSbr.^2.75 .*( abs(sind(2*WAVE.dPHIbr)).^0.5.*sign(WAVE.dPHIbr) - 2.*cosd(WAVE.dPHIbr).*dHS ); % added secondary component due to spatial wave height differences 
        QSmilhmass=TRANSP.qscal .* 0.15.*  TRANSP.rhos/(TRANSP.rhos-TRANSP.rhow) .* WAVE.TP.^0.89 .* TRANSP.tanbeta.^0.86 .* TRANSP.d50.^-0.69 .* WAVE.HSbr.^2.75 .*( abs(sind(2*WAVE.dPHIbr)).^0.5.*sign(WAVE.dPHIbr) - (2/TRANSP.tanbetasetup).*cosd(WAVE.dPHIbr).*dHS ); % added secondary component due to spatial wave height differences + the slope effect on this current
        QS = 365*24*60*60*(QSmilhmass / TRANSP.rhos) /(1.0-TRANSP.porosity);
        QS(abs(WAVE.dPHIbr)>90)=0;     % Set transport to 0 when angle exceeds 180 degrees
    
    %% Transport : Van Rijn (2014)        
    elseif ~isempty(findstr(lower(TRANSP.trform),'vr14'))
        vtide=0;
        kswell=0.015.*TRANSP.pswell+(1-0.01.*TRANSP.pswell);
        %vwave=0.3.*(TRANSP.g.*WAVE.HSbr).^0.5.*sind(2.*WAVE.dPHIbr); 
        %vwave=0.3.*(TRANSP.g.*WAVE.HSbr).^0.5.*(sind(2.*WAVE.dPHIbr)-2*cosd(WAVE.dPHIbr).* dHS); % added secondary component due to spatial wave height differences 
        vwave=0.3.*(TRANSP.g.*WAVE.HSbr).^0.5.*(sind(2.*WAVE.dPHIbr)-(2/TRANSP.tanbetasetup)*cosd(WAVE.dPHIbr).* dHS); % added secondary component due to spatial wave height differences         
        vwave(abs(WAVE.dPHIbr)>90)=0;
        vtotal=vwave+vtide;
        QSvr14mass=TRANSP.qscal .* 0.0006 .* kswell .* TRANSP.rhos .* TRANSP.tanbeta.^0.4 .*TRANSP.d50.^-0.6 .* WAVE.HSbr.^2.6 .* vtotal;
        QS = 365*24*60*60*(QSvr14mass ./ TRANSP.rhos) ./(1.0-TRANSP.porosity); 
        QS(abs(WAVE.dPHIbr)>90)=0;     % Set transport to 0 when angle exceeds 180 degrees
    
    elseif ~isempty(findstr(lower(TRANSP.trform),'ray'))
        QS = -WAVE.c1.*WAVE.dPHItdp .* exp(-(WAVE.c2.*WAVE.dPHItdp).^2) + WAVE.QSoffset;
        QS(abs(WAVE.dPHItdp)>90)=0;     % Set transport to 0 when angle exceeds 180 degrees
        dQSdPHI = WAVE.c1.*exp(-(WAVE.c2.*WAVE.dPHItdp).^2) .* ( -2.*(WAVE.c2.*WAVE.dPHItdp).^2 + 1. ) *180./pi; 

    elseif ~isempty(findstr(lower(TRANSP.trform),'tideprof'))
        alpha=1.0;    % Default dissipation parameter Baldock
        gamma=0.78;   % Default breaking parameter Baldock
        for i=1:nq
            Hrms0=WAVE.HStdp(1,i)/sqrt(2);
            Tp=WAVE.TP(1,i);
            theta0=WAVE.dPHItdp(1,i);
            Hrms0(abs(theta0)>90)=0;
            theta0(abs(theta0)>90)=0;
            [QS(1,i),Slong_mean,Slong,v,vt,vw,Hrms,h] = transport_tidewave(TIDE.eta(i,:),TIDE.detads(i,:),TIDE.phi(i,:),TIDE.ss(i), ...
                                                                         TIDE.Ttide,TIDE.nT,TIDE.k(i,:),TIDE.cf,TIDE.hmin, ...
                                                                         Hrms0,Tp,theta0,alpha,gamma, ...
                                                                         TRANSP.ks,TRANSP.d50,TRANSP.d90,TIDE.hclosure, ...
                                                                         TIDE.x,TIDE.zb,TRANSP.acal,TRANSP.n);
            if TRANSP.submerged==1
                [QSti,Dltri,zsi]=transport_groynesubmerged(Slong,vw,h,TIDE.zb,TIDE.x,TRANSP.aw,TRANSP.gamma,WAVE.HStdp(1,i));
                QSt(:,i)=QSti(:);
                Dltr(:,i)=Dltri(:);
                zs(:,i)=zsi(:);
                QS(1,i)=mean(QSti);
            end            
        end
    else
        fprintf(' Warning : No valid transport formulation defined!\n');
    end

    %% Export transport rates over cross-shore profile (if trform='TIDEWAVE')
    if TRANSP.submerged==1
        TRANSP.QSt=QSt;
        TRANSP.Dltr=Dltr;
        TRANSP.zs=zs;
    end

    %% Export transport
    TRANSP.(fieldQS)=QS;
    if strcmpi(fieldQS,'QS')
        TRANSP.debug.QS0=TRANSP.QS;
    end
end

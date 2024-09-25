function [GROYNE] = transport_bypass(TRANSP,WAVE,COAST,STRUC,GROYNE,TIDE)
% function [GROYNE] = transport_bypass(TRANSP,WAVE,COAST,STRUC,GROYNE,TIDE)
% 
% Computes the bypass rate at groynes if transport is towards the groyne.
%
% INPUT: 
%     TRANSP
%        .QS        : longshore transport rate [m3/yr]
%     WAVE
%        .HStdp     : Hs at toe of dynamic profile [m]
%     COAST
%        .x         : x coordinates of the active coastal section [m]
%        .y         : y coordinates of the active coastal section [m]
%     STRUC
%        .xhard     : x coordinates of all hard structures [m]
%        .yhard     : y coordinates of all hard structures [m]
%     GROYNE
%        .x         : x coordinates of groynes [m]
%        .y         : y coordinates of groynes [m]
%        .phigroyne : orientation of groyne [°N]
% 
% OUTPUT: 
%     GROYNE
%        .QS        : bypass transport at groynes [Mx2] array [m3/yr]
%        .tipindx   : index of most seaward point on each side of the structure [Mx2]
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

    if ~isempty(GROYNE.QS)

        %% determine cross-shore extent of groynes w.r.t. coastline
        for side=1:2
            idgroyne=find(GROYNE.idcoast(:,side)==COAST.i_mc); 
            GROYNE.QS(idgroyne,side)=0; % initialize at 0 
            
            if side==1
                % end point (left side of groyne & right side of coastal segment)       
                xbnd=COAST.x(end); 
                ybnd=COAST.y(end); 
                HSbnd=WAVE.HStdp(end); 
                if TRANSP.submerged==1
                    loca=size(TRANSP.QSt,1); % at end of coastal cell
                    QSbnd=TRANSP.QSt(loca,:);
                else
                    QSbnd=TRANSP.QS(end); 
                end
                iwave=COAST.nq;
            else
                % start point (right side of groyne & left side of coastal segment)
                xbnd=COAST.x(1); 
                ybnd=COAST.y(1); 
                HSbnd=WAVE.HStdp(1); 
                if TRANSP.submerged==1
                    loca=1; % at start of coastal cell 
                    QSbnd=TRANSP.QSt(loca,:);
                else
                    QSbnd=TRANSP.QS(1); 
                end
                iwave=1;
            end
            if (~isempty(idgroyne))
                ist=GROYNE.strucnum(idgroyne); 
                [x_str,y_str,n_str]=get_one_polygon(STRUC.xhard,STRUC.yhard,ist); 
                
                % furthest extent of idgroyne relative to shore normal 
                th=GROYNE.phigroyne(idgroyne);                          % orientation of groyne in °N
                proj=(x_str(2:end-1)-xbnd)*sind(th)+(y_str(2:end-1)-ybnd)*cosd(th); 
                [ystr,tipindx0]=max(proj);                              % ystr is the most seaward extent of the groyne on this side
                GROYNE.tipindx(idgroyne,side)=tipindx0+1;               % index of most seaward point on this side
                eps=1e-1; 
                ystr=max(ystr,eps); 
                
                % compute bypass fraction based on proxy for cross-shore extent of longshore transport and dean profile
                if (QSbnd>0 && side==1) || (QSbnd<0 && side==2)  
                    
                    Ap=(1.04+0.086*log(TRANSP.d50))^2;                    % sediment scale parameter
                    Ds=Ap*ystr^(2/3);Dsa=Ds;                              % depth at tip of structure according to dean profile
                    Dlt=TRANSP.aw.*WAVE.HStdp(iwave)/WAVE.gamma;          % most seaward extent of the longshore current
                    if TRANSP.submerged==1
                        % this is temporary unvalidated code for submerged groynes 
                        xx=interp1(zp,xp,0); % note that this only works if the profile is continiously going up in landward direction (i.e. without any humps / bars)
                        zz=interp1(xp,zp,xx-ystr);
                        if isnan(zz)
                            zz=min(zp);
                        end
                        Ds=-zz+TRANSP.zs(loca,:);
                        Dsa=TRANSP.zs(loca,:)-STRUC.groinelev(ist);
                        Dlt=TRANSP.Dltr(loca,:)+TRANSP.zs(loca,:);
                    end
                    
                    % option to use other proxies for the depth-of-closure
                    hin=2.28*WAVE.HStdp(iwave) - 68.5*(WAVE.HStdp(iwave).^2 ./(9.81.*WAVE.TP(iwave).^2)); % Inner closure-depth (Hallermeier, 1981)
                    hout=0.013*WAVE.HStdp(iwave)*WAVE.TP(iwave)*sqrt(9.81/(TRANSP.d50*1.65));
                    
                    % set the maxbypass ratio
                    % Bypass is restricted in case the coastline is further seaward at the downdrift side,
                    % The ratio of the effective wave angle over the structure is then used to scale the bypass.
                    % So, less than the maxbypass is used when other side is too much seaward (e.g. due to accretion against the leeward side)
                    maxbypass=max(TRANSP.bypasscontrfac,1);
                    if maxbypass>1
                        PHIcoast_groyne=mod(GROYNE.phicoast(idgroyne)-GROYNE.phigroyne(idgroyne)+180,360)-180;   % relative angle of the groyne orientation w.r.t. coast orientation (positive means transport towards the right)
                        PHI_wave_tdp=mod(WAVE.PHItdp(iwave)-GROYNE.phicoast(idgroyne)+180,360)-180;    % relative angle of the nearshore wave climate w.r.t. coast orientation (positive means transport towards the right)
                        effectivewaveangle = max( (PHI_wave_tdp - PHIcoast_groyne) / PHI_wave_tdp, 1);
                        maxbypass=min(maxbypass , effectivewaveangle);
                    end
                    
                    % Special case where transport is towards structure due to wave diffraction at the leeward side of the structure, 
                    % but in that case the bypass will not be enhanced. The offshore wave direction is used to identify this situation. 
                    PHI_groyne_o=mod(WAVE.PHIo(iwave)-GROYNE.phicoast(idgroyne)+180,360)-180;      % relative angle of the offshore wave climate w.r.t. coast orientation (positive means transport towards the right)
                    if (side==1 && PHI_groyne_o>0) || (side==2 && PHI_groyne_o<0)
                        maxbypass=max(1-abs(PHI_groyne_o),0); % decreases to 0 over the relative incoming angle of 1 to 0 degree, which is more stable than defining just 0 here
                    end
                    
                    % Determine the fraction of bypassing transport 'BPF' 
                    % using the ratio of the depth of the dean profile at the tip 'ds' and the depth-of-closure proxy 'dlt' of the longshore current.
                    BPF=maxbypass-Ds./Dlt;                               % fraction of bypassing sediment
                    
                    if TRANSP.submerged==1
                        DsaDs=min(max(Dsa./Ds,0),1);
                        redu=min(max(-0.9728.*min(DsaDs,0.8)+0.1296.*(Ds./Dlt)+1.0227),1);
                        id=find(DsaDs>0.8);
                        redu(id)=redu(id).*(1-DsaDs(id))/0.2;
                        BPF=maxbypass-Ds./Dlt.*redu;
                    end
                    if isnan(BPF)                                       % if aw is zero, switch off bypassing
                       BPF=0; 
                    end 
                    BPF=min(max(BPF,0),maxbypass); 
                    GROYNE.QS(idgroyne,side)=mean(BPF.*QSbnd);
                    % debug info
                    if 0
                        fprintf('%8.0f %8.0f %8.0f %8.3f %8.0f \n',idgroyne,ist,side,BPF,BPF*QSbnd)
                        if side==1
                        hold on;plot(x_str(GROYNE.tipindx(idgroyne,side)),y_str(GROYNE.tipindx(idgroyne,side)),'g*');
                        else
                        hold on;plot(x_str(GROYNE.tipindx(idgroyne,side)),y_str(GROYNE.tipindx(idgroyne,side)),'r*');
                        end
                    end
                end
            end
        end
    end
end


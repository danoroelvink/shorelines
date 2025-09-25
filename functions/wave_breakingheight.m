function [WAVE]=wave_breakingheight(WAVE,TRANSP)
% function [WAVE]=wave_breakingheight(WAVE,TRANSP)
% 
% This routine computes the refraction and shoaling from the nearshore point (TDP) to the point of breaking (BR).
% 
% INPUT: 
%    WAVE
%         .dnearshore        : depth at toe of dynamic profile [m]
%         .HStdp             : wave height at nearshore location (or diffracted wave height) [Nx1]
%         .TP                : wave period [s]
%         .dPHItdp           : relative angle of waves with respect to the coastline at nearshore location [°] 
%         .dPHIcrit          : Critical orientation of the coastline with respect to waves at the depth-of-closure (Nx2) [°]
%         .gamma             : depth-induced breaking coefficient [-]
%    TRANSP
%         .suppresshighangle : option to suppress high-angle instabilities by maximizing transport at angles beyond dphicrit (0/1)
%         .trform            : transport formulation (either 'CERC', 'KAMP', 'MILH', 'CERC3', 'VR14')
% 
% OUTPUT:
%    WAVE
%         .HSbr              : breaking wave height (or diffracted wave height) [Nx1]
%         .dPHIbr            : relative angle of waves with respect to the coastline at the point of breaking [°]
%         .dPHItdp           : relative angle of waves with respect to the coastline at nearshore location [°] (updated)
%         .hbr               : depth at point of breaking [m]
%         .cbr               : wave celerity at point of breaking [m/s]
%         .nbr               : deep/shallow water wave number [-]
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

    eps=1d-5;
    hmin=0.1;
    hbrmin=0.05; % minimum depthg of breaking
    WAVE.HSbr=zeros(size(WAVE.dPHItdp));
    WAVE.dPHIbr=zeros(size(WAVE.dPHItdp));
    WAVE.hbr=zeros(size(WAVE.dPHItdp));
    WAVE.cbr=zeros(size(WAVE.dPHItdp));
    WAVE.nbr=zeros(size(WAVE.dPHItdp));
    
    %% Limit high-angle instabilities by forcing maximum transport when the angle is larger than dPHIcrit
    if TRANSP.suppresshighangle==1 
        if isfield(WAVE,'dPHIcrit_mc0')
            dPHIcrit = get_one_polygon( WAVE.dPHIcrit_mc0,WAVE.i_mc);
            if length(dPHIcrit(~isnan(dPHIcrit)))==length(WAVE.dPHItdp)
                WAVE.dPHItdp=sign(WAVE.dPHItdp).*min([abs(WAVE.dPHItdp);dPHIcrit],[],1);        
            elseif length(dPHIcrit(~isnan(dPHIcrit)))>=1
                WAVE.dPHItdp=sign(WAVE.dPHItdp).*min(abs(WAVE.dPHItdp),median(dPHIcrit));
            else
                WAVE.dPHItdp=sign(WAVE.dPHItdp).*min(abs(WAVE.dPHItdp),45);    
            end
        else
            WAVE.dPHItdp=sign(WAVE.dPHItdp).*min(abs(WAVE.dPHItdp),45);        
        end
    end
    
    if ~strcmpi(TRANSP.trform,'RAY') && ~strcmpi(TRANSP.trform,'CERC') && ~strcmpi(TRANSP.trform,'CERC2')
        for i=1:size(WAVE.HStdp,1)   
        for j=1:size(WAVE.HStdp,2)
            [~,ctdp,~,ntdp]=get_disper(WAVE.dnearshore,WAVE.TP(i,j));
            cosPHI=max(cosd(WAVE.dPHItdp(i,j)),eps);
            sinPHI=sind(WAVE.dPHItdp(i,j));
            iter=1; % First estimate: Hbr=WAVE.HStdp
            hbr1=max(WAVE.HStdp(i,j)./WAVE.gamma,eps);
            [hbr2,dPHIbr0,err1,cbr0,nbr0]=wave_shoalref(hbr1,WAVE.TP(i,j),WAVE.gamma,WAVE.HStdp(i,j),ctdp,ntdp,sinPHI,cosPHI);
            hbr2b=max(hbr2,eps);
            iter=2; % Second estimate: fill in hbr2
            [hbr3,dPHIbr0,err2,cbr0,nbr0]=wave_shoalref(hbr2b,WAVE.TP(i,j),WAVE.gamma,WAVE.HStdp(i,j),ctdp,ntdp,sinPHI,cosPHI);
            
            if abs(err2)<eps
                hbrnew=hbr3;
                errnew=err2;
            else
                for iter=3:10 % Following estimates: inter/extrapolate from last two WAVE.hbr/err pairs
                    hbrnew=max(hbr1-err1*(hbr2-hbr1)/(err2-err1),eps);
                    [hbrest,dPHIbr0,errnew,cbr0,nbr0]=wave_shoalref(hbrnew,WAVE.TP(i,j),WAVE.gamma,WAVE.HStdp(i,j),ctdp,ntdp,sinPHI,cosPHI);
                    if abs(errnew)>eps
                        hbr1=hbr2;err1=err2;
                        hbr2=hbrnew;err2=errnew;
                    else
                        break
                    end
                end
            end
            
            %disp([num2str(PHI(iphi),'%4.1f'),' ',num2str(iter,'%4i'),' ',num2str(hbrnew,'%4.2f'),' ',num2str(errnew,'%6.4f')])
            hbrnew(isnan(hbrnew))=0;
            WAVE.HSbr(i,j)=hbrnew*WAVE.gamma;
            WAVE.dPHIbr(i,j)=dPHIbr0;
            WAVE.dPHIbr(isnan(WAVE.dPHIbr))=0;
            WAVE.hbr(i,j)=hbrnew;
            WAVE.cbr(i,j)=cbr0;
            WAVE.nbr(i,j)=nbr0;
        end
        end
           
    else 
        % in case TRANSP.trform is 'RAY', 'CERC' or 'CERC2'
        WAVE.HSbr=WAVE.HStdp;
        WAVE.dPHIbr=WAVE.dPHItdp;
        WAVE.hbr=max(WAVE.HStdp./WAVE.gamma,0.1*hbrmin);
        [~,cbr,~,nbr]=get_disper(WAVE.hbr,WAVE.TP);
        WAVE.cbr=cbr;
        WAVE.nbr=nbr;
    end
    
end

function [QSt,Dltr,zs]=transport_groynesubmerged(Slong,vw,h,zb,x,Aw,gammabr,HStdp)
% function [QSt,Dltr,zs]=transport_groynesubmerged(Slong,vw,h,zb,x,Aw,gammabr,HStdp)
% 
% Transport_groynesubmerged computes the tide-wave driven transports over a submerged groyne.
% It should be used in combination with the TIDEPROF formulation. 
% 
% INPUT: 
%     Slong       : longshore transport over the cross-shore profile for all tide timesteps [nx nt] [m3/m/yr]
%     vw          : longshore wave-driven current velocity, 2D matrix [nx nT] [m/s]
%     h           : depth over the cross-shore profile [nx nT] [m]
%     zb          : z-coordinates of cross-shore profiles [m]
%     x           : x-coordinates of cross-shore profiles [m]
%     Aw          : coefficient for cross-shore distribution of alongshore transport at groynes [-]
%     gammabr     : depth-induced breaking coefficient [-]
%     HStdp       : Wave height at nearshore location (or diffracted wave height) [m]
% 
% OUTPUT:
%     QSt         : adjusted total longshore transport over submerged groyne [m3/yr]
%     Dltr        : critical depth at tip of groyne for bypassing, which is corrected here [m]
%     zs          : water depth at submerged groyne
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

    for tt=1:24
        % if tide is dominant, the max can be at the offshore end. 
        % Transport is in that case anyway small, so not that important
        noDLT=0;
        if max(vw(:,tt))>0
            slongmax=max((Slong(:,tt)));
        else
            slongmax=abs(min((Slong(:,tt))));
        end
        [ind]=find(abs(Slong(:,tt))==slongmax);
        jj=ind(1);
        slongmaxe=Slong(jj,tt);
    
        while abs(Slong(jj,tt))>0.1*slongmax & (noDLT==0 & slongmaxe*Slong(jj,tt)>0)
            jj=jj-1;                         
            if jj<3
                noDLT=1;
                jj=1;
            end
        end
        if noDLT==0
            Dltr(tt)=-zb(jj);
        else
            Dltr(tt)=Aw*HStdp/gammabr;
        end
        QSt(tt)=sum(Slong(:,tt))*(x(2)-x(1));
        zs(tt)=h(1,tt)+zb(1);
    end
end


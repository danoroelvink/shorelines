function [hbrnew,dPHIbr,err,cbr,nbr]=wave_shoalref(hbr,tper,gamma,hstdp,ctdp,ntdp,sinPHIw,cosPHIw)
% function [hbrnew,dPHIbr,err,cbr,nbr]=wave_shoalref(hbr,tper,gamma,hstdp,ctdp,ntdp,sinPHIw,cosPHIw)
%
% This routine computes the change in wave height and breaking depth as a 
% result of shoaling of the waves and wave-height reduction due to refraction.
% 
% INPUT: 
%     hbr               : estimate of depth at point of breaking [m]
%     tper              : wave period [s]
%     gamma             : wave breaking parameter gamma [-]
%     hstdp             : significant wave height at depth-of-closure [m]
%     ctdp              : wave celerity at depth-of-closure [m/s]
%     ntdp              : wave number at depth-of-closure [-]
%     sinPHIw           : sinus component of incoming wave angle
%     cosPHIw           : cosine component of incoming wave angle
% 
% OUTPUT:
%     hbrnew            : depth at point of breaking [m] (new iteration)
%     dPHIbr            : relative wave angle w.r.t. coast at point of breaking [°]
%     err               : relative change in depth at point of breaking w.r.t. previous iteration (hbrnew-hbr) [m]
%     cbr               : wave celerity at point of breaking [m/s]
%     nbr               : deep/shallow water wave number [-]
%
% Copyright notice
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

    [~,cbr,~,nbr]=get_disper(hbr,tper);
    if abs(cbr./ctdp.*sinPHIw)<1
        dPHIbr=asind(cbr./ctdp.*sinPHIw);
        if (nbr.*cbr.*cosd(dPHIbr))>0 && abs(cosd(dPHIbr))>1e-3
            hstdpbr=hstdp.*sqrt(ntdp.*ctdp.*cosPHIw./(nbr.*cbr.*cosd(dPHIbr)));
            %hstdpbr=hstdp;
        else
            dPHIbr=acosd(cosPHIw); % take the input wave angle as representative for the relative incidence angle at breaking
            hstdpbr=hstdp;  % not so relevant case where the wave approaches very obliquely, or depth is 0, likely very small TP of hstdp
        end
    else
        dPHIbr=acosd(cosPHIw); % take the input wave angle as representative for the relative incidence angle at breaking
        hstdpbr=hstdp; % in case we have a extremely oblique wave incidence (close to 90°) in combination with further offshore depth of hbr estimate (compared to old 'hbr'), likely very small TP of hstdp
    end

    hbrnew=hstdpbr/gamma;
    err=hbrnew-hbr;
end

function [hbrnew,dPHIbr,err,cbr,nbr]=wave_shoalref(hbr,tper,gamma,HS,ctdp,ntdp,sinPHIw,cosPHIw)
% function [hbrnew,dPHIbr,err,cbr,nbr]=wave_shoalref(hbr,tper,gamma,HS,ctdp,ntdp,sinPHIw,cosPHIw)
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

    [khb,cbr]=wave_GUO2002(tper,hbr);
    nbr = 0.5.*(1.+(2.*khb)./sinh(2.*khb));
    if abs(cbr./ctdp.*sinPHIw)<1
        dPHIbr=asind(cbr./ctdp.*sinPHIw);
        if (nbr.*cbr.*cosd(dPHIbr))>0 && abs(cosd(dPHIbr))>1e-3
            HSbr=HS.*sqrt(ntdp.*ctdp.*cosPHIw./(nbr.*cbr.*cosd(dPHIbr)));
            %HSbr=HS;
        else
            dPHIbr=acosd(cosPHIw); % take the input wave angle as representative for the relative incidence angle at breaking
            HSbr=HS;  % not so relevant case where the wave approaches very obliquely, or depth is 0, likely very small TP of HS
        end
    else
        dPHIbr=acosd(cosPHIw); % take the input wave angle as representative for the relative incidence angle at breaking
        HSbr=HS; % in case we have a extremely oblique wave incidence (close to 90°) in combination with further offshore depth of hbr estimate (compared to old 'hbr'), likely very small TP of HS
    end

    hbrnew=HSbr/gamma;
    err=hbrnew-hbr;
end

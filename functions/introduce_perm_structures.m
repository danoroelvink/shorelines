function [WAVE]=introduce_perm_structures(STRUC,WAVE,COAST)
% [Hsout]=introduce_perm_structures(STRUC,WAVE,COAST)
% 
% INPUT:
%   STRUCT
%       .x_perm         Coordinates of permeable breakwaters     
%       .y_perm         Coordinates of permeable breakwaters     
%   WAVE
%       .HStdp          Significant wave height at considered time instance (m)
%       .wavetransm     Array with transmission coefficients per section
%   COAST
%       .x              Coordinates of current coast section
%       .y              Coordinates of current coast section
% 
% OUTPUT:
%   WAVE.HSo            Transmitted significant wave height at considered time instance at offshore location (m)
%   WAVE.HStdp          Transmitted significant wave height at considered time instance at depth-of-closure (m)
%   WAVE.HSbr           Transmitted significant wave height at considered time instance at point of breaking (m)
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

    if STRUC.perm
        [x_p,y_p,n_mc,i1,i2]=get_one_polygon(STRUC.x_perm,STRUC.y_perm,1);
        for i_mc=1:n_mc
            [x_p,y_p,n_mc,i1,i2]=get_one_polygon(STRUC.x_perm,STRUC.y_perm,i_mc);
            [ shadowP ] = find_shadows_mc(COAST.x,COAST.y,x_p,y_p,WAVE.PHIo,WAVE.PHItdp,WAVE.PHIbr);
            WAVE.HSo(shadowP)=WAVE.HSo(shadowP)*WAVE.wavetransm(i_mc);
            WAVE.HStdp(shadowP)=WAVE.HStdp(shadowP)*WAVE.wavetransm(i_mc);
            WAVE.HSbr(shadowP)=WAVE.HSbr(shadowP)*WAVE.wavetransm(i_mc);
        end
    end
end

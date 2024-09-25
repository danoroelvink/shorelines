function [WAVE]=introduce_perm_structures(STRUC,WAVE,COAST)
% function [WAVE]=introduce_perm_structures(STRUC,WAVE,COAST)
% 
% This function introduces permeable structures by adjusting the 
% wave conditions at the rear side of the structure. 
%
% INPUT: 
%   STRUCT
%       .xperm       : x-coordinates of permeable breakwaters     
%       .yperm       : y-coordinates of permeable breakwaters     
%   WAVE
%       .HSo         : significant wave height at offshore location at considered time instance (m)
%       .HStdp       : significant wave height at the depth-of-closure at considered time instance (m)
%       .HSbr        : significant wave height at the point of breaking at considered time instance (m)
%       .PHIo        : wave direction at offshore location at considered time instance (°N)
%       .PHItdp      : wave direction at the depth-of-closure at considered time instance (°N)
%       .PHIbr       : wave direction at the point of breaking at considered time instance (°N)
%       .wavetransm  : Array with transmission coefficients per section
%   COAST
%       .x           : x-coordinates of current coast section
%       .y           : y-coordinates of current coast section
% 
% OUTPUT:
%   WAVE.HSo         : transmitted significant wave height at considered time instance at offshore location (m)
%   WAVE.HStdp       : transmitted significant wave height at considered time instance at depth-of-closure (m)
%   WAVE.HSbr        : transmitted significant wave height at considered time instance at point of breaking (m)
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

    if STRUC.perm
        [x_p,y_p,n_mc,i1,i2]=get_one_polygon(STRUC.xperm,STRUC.yperm,1);
        for i_mc=1:n_mc
            [x_p,y_p,n_mc,i1,i2]=get_one_polygon(STRUC.xperm,STRUC.yperm,i_mc);
            [ shadowP ] = find_shadows_mc(COAST.x,COAST.y,x_p,y_p,WAVE.PHIo,WAVE.PHItdp,WAVE.PHIbr);
            WAVE.HSo(shadowP)=WAVE.HSo(shadowP)*WAVE.wavetransm(i_mc);
            WAVE.HStdp(shadowP)=WAVE.HStdp(shadowP)*WAVE.wavetransm(i_mc);
            WAVE.HSbr(shadowP)=WAVE.HSbr(shadowP)*WAVE.wavetransm(i_mc);
        end
    end
end

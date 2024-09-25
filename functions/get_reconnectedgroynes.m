function [COAST]=get_reconnectedgroynes(COAST,GROYNE)   
% function [COAST]=get_reconnectedgroynes(COAST,GROYNE)   
%
% INPUT:  
%    COAST       : coast data strcuture (e.g. field 'x_mc')
%    WAVE        : wave data structure (e.g. field 'Hso_mc')
%    TRANSP      : transport properties (e.g. field 'QS_mc')
%    GROYNE      : groyne data structure with field 'idcoast'   
%
% OUTPUT:
%    COAST       : coast data strcuture (e.g. field 'x_mc')
%    WAVE        : wave data structure (e.g. field 'Hso_mc')
%    TRANSP      : transport properties (e.g. field 'QS_mc')
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2020 Deltares & IHE-Delft
%
%       Bas Huisman
%       bas.huisman@deltares.nl
%       Boussinesqweg 1
%       2629HV Delft
%
%       Dano Roelvink
%       d.roelvink@un-ihe.org
%       Westvest 7
%       2611AX Delft
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

    COAST0=COAST;
    if isfield(GROYNE,'idcoast')
        for ig=1:size(GROYNE.idcoast,1)
            % get the segment indices at both sides (from GROYNE.idcoast)
            % obtain the COAST, WAVE and TRANSP information of these segments
            i_mc1=GROYNE.idcoast(ig,1);
            i_mc2=GROYNE.idcoast(ig,2);
            if ~isnan(i_mc1)
            [x1,n_mc1,i1start,i1end ] = get_one_polygon( COAST0.x_mc,i_mc1);
            end
            if ~isnan(i_mc2)
            [x2,n_mc2,i2start,i2end ] = get_one_polygon( COAST0.x_mc,i_mc2);
            end
            
            % join the 2 segments at both sides of the groyne with an additional point in-between
            if isnan(i_mc1) || isnan(i_mc2)
                % one-sided groyne at edge of model
            elseif i2start-1==i1end+1
                % groyne somewhere inside the model
                igroyne = i1end+1;
                %xg=0.5*(COAST.x_mc(i1end)+COAST.x_mc(i2start)); % old method of point exactly in-between the two points.
                %yg=0.5*(COAST.y_mc(i1end)+COAST.y_mc(i2start));
                % new method which uses the stored first point that is inside the groyne that is found from the first time the coast intersected the structure.
                % this avoids problems with points at the boundary of the groyne. 
                xg=GROYNE.xg(ig); 
                yg=GROYNE.yg(ig); 
                
                COAST.x_mc=[COAST.x_mc(1:i1end),xg,COAST.x_mc(i2start:end)]; 
                COAST.y_mc=[COAST.y_mc(1:i1end),yg,COAST.y_mc(i2start:end)]; 
            else
                disp('Can not re-connect groyne beacuse the coastal segments are not adjacent!')
            end
        end
        %% make transport points xq_mc and yq_mc
        [COAST]=get_transportpoints(COAST);
    end
end
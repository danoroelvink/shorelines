function [COAST]=get_coastline_orientation(COAST)
% function [COAST]=get_coastline_orientation(COAST)
% 
% INPUT:
%   COAST
%     .x          x-coordinate of coastal segment [m]
%     .y          y-coordinate of coastal segment [m]
%
% OUTPUT:
%   COAST
%     .PHIc       Coastline orientation at transport grid points [°N]
%     .PHIcxy      Coastline orientation at coastline grid points [°N]
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2021 IHE Delft & Deltares
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

    % coast angles (PHIc) at transport points (only in-between coastline points, not at 'boundary transport points')    
    % angle in-between grid cell points
    COAST.PHIc=mod(360.0-atan2d(diff(COAST.y),diff(COAST.x)),360);
    
    % coast angles (PHIcxy) at coastline points
    COAST.PHIcxy=mod(360.0-atan2d(diff(COAST.yq),diff(COAST.xq)),360);
    
    % store initial coastline
    if isnan(COAST.PHIc0bnd(1))
        COAST.PHIc0bnd(1)=COAST.PHIc(1);
    end
    if isnan(COAST.PHIc0bnd(2))
        COAST.PHIc0bnd(2)=COAST.PHIc(end);
    end
       
    % add extra boundary condition cell    
    if COAST.cyclic
        % make sure to use end and start points to make it cyclical
        COAST.PHIc=[COAST.PHIc(end),COAST.PHIc,COAST.PHIc(1)];
    else
        % extend coast-angles in case of an open section
        COAST.PHIc=[COAST.PHIc(1),COAST.PHIc,COAST.PHIc(end)];
        
        % store coast-angle at bounary at t0 for 'angleconstant' boundary conditions
        % or use a prescribed angle for the boundary
        % and re-use this angle at later timesteps. 
        if strcmpi(COAST.boundary_condition_start{1},'Angleconstant') || ...
           strcmpi(COAST.boundary_condition_start{1},'Gradient') 
            if isnan(COAST.boundary_condition_start{2})
                COAST.boundary_condition_start{2}=COAST.PHIc(1);
            else
                COAST.PHIc(1)=COAST.boundary_condition_start{2};
            end
        end
        if strcmpi(COAST.boundary_condition_end{1},'Angleconstant') || ...
           strcmpi(COAST.boundary_condition_end{1},'Gradient') 
            if isnan(COAST.boundary_condition_end{2})
                COAST.boundary_condition_end{2}=COAST.PHIc(end);
            else
                COAST.PHIc(end)=COAST.boundary_condition_end{2};
            end
        end
    end
    
    % smooth coastline angles
    % the smoothed coastline PHIcs is only used for the refraction in the nearshore, 
    % and not for computing the transport rates. 
    COAST.PHIcs=COAST.PHIc;
    if COAST.smoothrefrac>0 
        COAST.PHIcs=get_smoothdata(COAST.PHIc,'Angle',min(COAST.smoothrefrac,1)); 
    end
end
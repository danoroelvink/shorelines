function [COAST]=get_coastline_orientation(COAST)
% function [COAST]=get_coastline_orientation(COAST)
% 
% INPUT: 
%   COAST
%     .x                      : x-coordinates at coastline points [m]
%     .y                      : y-coordinates at coastline points [m]
%     .xq                     : x-coordinates at qs-points [m]
%     .yq                     : y-coordinates at qs-points [m]
%     .cyclic                 : index indicating cyclicality of the active element (0/1)
%
% OUTPUT:
%   COAST
%     .PHIc                   : coastline orientation at transport grid points [°N]
%     .PHIcxy                 : coastline orientation at coastline grid points [°N]
%     .PHIcs                  : smoothed coastline orientation at transport grid points [°N], used only for refraction
%     .PHIc0bnd               : orientation of the coastline at the boundary points at t0 [°N]
%     .boundaryconditionstart : orientation is fixed of the boundary condition at to when 'Angleconstant'
%     .boundaryconditionend   : orientation is fixed of the boundary condition at to when 'Angleconstant'
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
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses>
%   --------------------------------------------------------------------

    % coast angles (PHIc) at transport points  
    % angle in-between grid cell points (not at 'boundary transport points')   
    COAST.PHIc=mod(360.0-atan2d(diff(COAST.y),diff(COAST.x)),360);

    % coast angles (PHIcxy) at coastline points
    COAST.PHIcxy=mod(360.0-atan2d(diff(COAST.yq),diff(COAST.xq)),360);

    % add extra boundary condition cell    
    if COAST.cyclic
        % coast angles (PHIc) at transport points 
        COAST.PHIc=[COAST.PHIc(end),COAST.PHIc,COAST.PHIc(1)];
        
    else
        % extend coast-angles in case of an open section
        COAST.PHIc=[COAST.PHIcxy(1),COAST.PHIc,COAST.PHIcxy(end)];
        
        % store coast-angle at bounary at t0 for 'angleconstant' boundary conditions
        % or use a prescribed angle for the boundary
        % and re-use this angle at later timesteps. 
        if strcmpi(COAST.boundaryconditionstart{1},'Angleconstant') || ...
           strcmpi(COAST.boundaryconditionstart{1},'Gradient') 
            if isnan(COAST.boundaryconditionstart{2})
                COAST.boundaryconditionstart{2}=COAST.PHIc(1);
            else
                COAST.PHIc(1)=COAST.boundaryconditionstart{2};
            end
        end
        if strcmpi(COAST.boundaryconditionend{1},'Angleconstant') || ...
           strcmpi(COAST.boundaryconditionend{1},'Gradient') 
            if isnan(COAST.boundaryconditionend{2})
                COAST.boundaryconditionend{2}=COAST.PHIc(end);
            else
                COAST.PHIc(end)=COAST.boundaryconditionend{2};
            end
        end
    end
    
    % store initial coastline
    if isnan(COAST.PHIc0bnd(1))
        COAST.PHIc0bnd(1)=COAST.PHIc(1);
    end
    if isnan(COAST.PHIc0bnd(2))
        COAST.PHIc0bnd(2)=COAST.PHIc(end);
    end
    if COAST.gridchange==1 || COAST.mergegrid==1
        COAST.PHIcxy0=COAST.PHIcxy;
        if isfield(COAST,'PHIcxy1_mc')
        method='weighted_distance';
        try
        [~,COAST.PHIcxy0]=get_interpolation_on_grid(method,COAST.x,COAST.y,COAST.x1_mc,COAST.y1_mc,'',COAST.PHIcxy1_mc);
        end
        end
    else
        COAST.PHIcxy0=get_one_polygon(COAST.PHIcxy1_mc,COAST.i_mc);
    end
    
    % smooth coastline angles
    % the smoothed coastline PHIcs is only used for the refraction in the nearshore, 
    % and not for computing the transport rates. 
    COAST.PHIcs=COAST.PHIc;
    if COAST.smoothrefrac>0 
        COAST.PHIcs=get_smoothdata(COAST.PHIc,'Angle',min(COAST.smoothrefrac,1)); 
    end
end
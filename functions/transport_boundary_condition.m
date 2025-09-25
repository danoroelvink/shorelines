function [TRANSP]=transport_boundary_condition(TRANSP,COAST,GROYNE)
% function [TRANSP]=transport_boundary_condition(TRANSP,COAST,GROYNE)
% 
% This function applies boundary conditions to non cyclic shorelines. 
% Cyclical elements are dealt with in the coastline_change function. 
% Type of boundaries are:
%   'Fixed' or 'neumann'           Neumann boundary with fixed coastline position
%   'Periodic'                     Periodic boundary with synchronized Q + position change at both ends of the coastal element
%   'Closed'                       Dirichlet or 'wall boundary' (transport can be enforced using {'Closed',QSrate} )
%   'Angleconstant' or 'Gradient'  Constant orientation/angle of the boundary points (coast angle can be enforced using {'Angleconstant',Angle} )
%                                  The angle is prescribed in get_coastline_orientation for the last QS index
% 
% INPUT: 
%     TRANSP
%        .QS                     : longshore transport rate [m3/yr]
%        .boundaryconditionstart : boundary condition at start of open coastline element (see types in description, e.g. 'Fixed' or {'Closed',5000} )
%        .boundaryconditionend   : boundary condition at end of open coastline element (see types in description, e.g. 'Fixed' or {'Closed',5000} )
%     WAVE
%        .HStdp                  : Hs at toe of dynamic profile [m]
%     COAST
%        .i_mc                   : index of active coastal element
%        .cyclic                 : index with cyclicality of the coastal element (0/1)
%        .BNDgroyne              : index showing the presence of structures at the start-end of each coastal section, and their structure id [number of coastal elements x 2-sides]
%     GROYNE
%        .idcoast                : index of segments at both sides of the groyne, used for recombining elements after the coastline change [Mx2]
% 
% OUTPUT: 
%     TRANSP
%        .QS                     : longshore transport rate, updated at boundaries [m3/yr]
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
    
    %% ADD BOUNDARY CONDITIONS
    QS0=TRANSP.QS;
    if ~COAST.cyclic
        if COAST.BNDgroyne(COAST.i_mc,1)~=0
            % transport has been computed in the bypassing routine 
            side=2;
            idgroyne=find(GROYNE.idcoast(:,side)==COAST.i_mc); 
            if (~isempty(idgroyne))
                %TRANSP.QS(:,1)=min(TRANSP.QS(:,1),GROYNE.QS(idgroyne,side));
            end
        else
            if strcmpi(TRANSP.boundaryconditionstart{1},'Periodic') && ~COAST.cyclic
                % periodic b.c.
                TRANSP.QS(:,1)=(QS0(:,end-1)+QS0(:,2))/2;
            end
            if strcmpi(TRANSP.boundaryconditionstart{1},'Closed') && ~COAST.cyclic
                % specify transport b.c.
                if isnan(TRANSP.boundaryconditionstart{2})
                    TRANSP.QS(:,1)=0;
                else
                    TRANSP.QS(:,1)=TRANSP.boundaryconditionstart{2};
                end
            end
            if (strcmpi(TRANSP.boundaryconditionstart{1},'Neumann') || strcmpi(TRANSP.boundaryconditionstart{1},'Fixed')) && ~COAST.cyclic
                % fixed coastline position b.c.
                TRANSP.QS(:,1)=TRANSP.QS(:,2);
            end
            if (strcmpi(TRANSP.boundaryconditionstart{1},'Angleconstant') || strcmpi(TRANSP.boundaryconditionstart{1},'Gradient')) && ~COAST.cyclic 
                % coast-angle constant b.c.
                % the angle is prescribed in get_coastline_orientation for the first QS index
            end
        end        
        
        if COAST.BNDgroyne(COAST.i_mc,2)~=0
            % transport has been computed in the bypassing routine 
            side=1;
            idgroyne=find(GROYNE.idcoast(:,side)==COAST.i_mc); 
            if (~isempty(idgroyne))                
                %TRANSP.QS(:,end)=max(TRANSP.QS(:,end),GROYNE.QS(idgroyne,side));
            end
        else
            if strcmpi(TRANSP.boundaryconditionend{1},'Periodic') && ~COAST.cyclic
                % periodic b.c.
                TRANSP.QS(:,end)=(QS0(:,end-1)+QS0(:,2))/2;
            end
            if strcmpi(TRANSP.boundaryconditionend{1},'Closed') && ~COAST.cyclic
                % specify transport b.c.
                if isnan(TRANSP.boundaryconditionend{2})
                    TRANSP.QS(:,end)=0;
                else
                    TRANSP.QS(:,end)=TRANSP.boundaryconditionend{2};
                end
            end
            if (strcmpi(TRANSP.boundaryconditionend{1},'Neumann') || strcmpi(TRANSP.boundaryconditionend{1},'Fixed')) && ~COAST.cyclic
                % fixed coastline position b.c.
                TRANSP.QS(:,end)=TRANSP.QS(:,end-1);
            end
            if (strcmpi(TRANSP.boundaryconditionend{1},'Angleconstant') || strcmpi(TRANSP.boundaryconditionend{1},'Gradient')) && ~COAST.cyclic 
                % coast-angle constant b.c.
                % The angle is prescribed in get_coastline_orientation for the last QS index
            end
        end
    end
end

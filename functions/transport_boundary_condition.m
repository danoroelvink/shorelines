function [TRANSP]=transport_boundary_condition(TRANSP,COAST,GROYNE)
% this function applied for non cyclic shorelines 
% In case of more than one (non-cyclic) section with different...
% ... boundary the code should be adjusted 
%
%'neumann'...Neumann boundary
%'CTAN'...Constant orientation 
%'func'...Dirichlet boundary 
%'closed'...Dirichlet (wall boundary)
%'periodic'...Periodic boundary Q+position   % testing phase 
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


    
    %% ADD BOUNDARY CONDITIONS
    QS0=TRANSP.QS;
    if ~COAST.cyclic
        if COAST.BNDgroyne(COAST.i_mc,1)~=0
            % transport has been computed in the bypassing routine 
            side=2;
            idgroyne=find(GROYNE.idcoast(:,side)==COAST.i_mc); 
            if (~isempty(idgroyne))
            %TRANSP.QS(1)=min(TRANSP.QS(1),GROYNE.QS(idgroyne,side));
            end
        else
            if strcmpi(TRANSP.boundary_condition_start{1},'Periodic') && ~COAST.cyclic
                % periodic b.c.
                TRANSP.QS(1)=(QS0(end-1)+QS0(2))/2;
            end
            if strcmpi(TRANSP.boundary_condition_start{1},'Closed') && ~COAST.cyclic
                % specify transport b.c.
                if isnan(TRANSP.boundary_condition_start{2})
                    TRANSP.QS(1)=0;
                else
                    TRANSP.QS(1)=TRANSP.boundary_condition_start{2};
                end
            end
            if (strcmpi(TRANSP.boundary_condition_start{1},'Neumann') || strcmpi(TRANSP.boundary_condition_start{1},'Fixed')) && ~COAST.cyclic
                % fixed coastline position b.c.
                TRANSP.QS(1)=TRANSP.QS(2);
            end
            if (strcmpi(TRANSP.boundary_condition_start{1},'Angleconstant') || strcmpi(TRANSP.boundary_condition_start{1},'Gradient')) && ~COAST.cyclic 
                % coast-angle constant b.c.
                % the angle is prescribed in get_coastline_orientation for the first QS index
            end
        end        
        
        if COAST.BNDgroyne(COAST.i_mc,2)~=0
            % transport has been computed in the bypassing routine 
            side=1;
            idgroyne=find(GROYNE.idcoast(:,side)==COAST.i_mc); 
            if (~isempty(idgroyne))                
            %TRANSP.QS(end)=max(TRANSP.QS(end),GROYNE.QS(idgroyne,side));
            end
        else
            if strcmpi(TRANSP.boundary_condition_end{1},'Periodic') && ~COAST.cyclic
                % periodic b.c.
                TRANSP.QS(end)=(QS0(end-1)+QS0(2))/2;
            end
            if strcmpi(TRANSP.boundary_condition_end{1},'Closed') && ~COAST.cyclic
                % specify transport b.c.
                if isnan(TRANSP.boundary_condition_end{2})
                    TRANSP.QS(end)=0;
                else
                    TRANSP.QS(end)=TRANSP.boundary_condition_end{2};
                end
            end
            if (strcmpi(TRANSP.boundary_condition_end{1},'Neumann') || strcmpi(TRANSP.boundary_condition_end{1},'Fixed')) && ~COAST.cyclic
                % fixed coastline position b.c.
                TRANSP.QS(end)=TRANSP.QS(end-1);
            end
            if (strcmpi(TRANSP.boundary_condition_end{1},'Angleconstant') || strcmpi(TRANSP.boundary_condition_end{1},'Gradient')) && ~COAST.cyclic 
                % coast-angle constant b.c.
                % the angle is prescribed in get_coastline_orientation for the last QS index
            end
        end
    end

%     % set nans to 0
%     if find(isnan(TRANSP.QS))>0
%         TRANSP.QS(find(isnan(TRANSP.QS)))=0;
%     end

end


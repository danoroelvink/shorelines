function [x_perm,y_perm]=prepare_permstructures(S)
% function [x_perm,y_perm]=prepare_perm_structures(S)
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

    if S.perm
        if isfield(S,'x_perm')
            x_perm=S.x_perm;
            y_perm=S.y_perm;
        elseif ~isempty(S.LDBstructures)
            xy_perm=load(S.LDBstructures);
            x_perm=xy_perm(:,1)'-S.XYoffset(1);
            y_perm=xy_perm(:,2)'-S.XYoffset(2);
            S.x_perm=x_perm;
            S.y_perm=y_perm;
        else
            figure(11);
            axis equal;
            xl=xlim;yl=ylim;
            htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add permeable structure (LMB); Next structure (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
            [x_perm,y_perm]=select_multi_polygon('k');
            set(htxt2,'Visible','off');
            S.x_perm=x_perm;
            S.y_perm=y_perm;
        end
    else
        x_perm=[];
        y_perm=[];
    end
    
end    

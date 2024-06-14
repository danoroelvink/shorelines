function [x_trans,y_trans,n_trans]=prepare_transects(S,x_mc,y_mc,x_hard,y_hard,x_dune,y_dune,xl,yl)
% function [x_trans,y_trans,n_trans]=prepare_transects(S,x_mc,y_mc,x_hard,y_hard,x_dune,y_dune,xl,yl)
% 
% a function to plot transects of cross-shore distances for both coastline and dune foot to a reference line...
% at specific times and transecs
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

    if  ~isempty(S.transect)
        xy_trans=load(S.transect);
        x_trans=xy_trans(:,1);
        y_trans=xy_trans(:,2);
        n_trans=length(x_trans)/2;
    else
        figure(12);
        plot(x_mc,y_mc,'k');hold on;
        if S.struct
            plot(x_hard,y_hard,'k');
        end
        if S.dune
            plot(x_dune,y_dune,'b');
            plot(xl,yl,'r');
            legend('Coastline','Structures','Dune Foot','Reference Line')
        end
        if S.dune && S.struct
            legend('Coastline','Structures','Dune Foot','Reference Line')
        elseif S.dune && ~S.struct
            legend('Coastline','Dune Foot','Reference Line')
        end
        
        axis equal;
        xlim(S.xlimits);
        ylim(S.ylimits);
        xlm=xlim;ylm=ylim;
        htxt2=text(xlm(1)+0.02*diff(xlm),ylm(2)-0.01*diff(ylm),'Add transect (LMB); Next transect (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.75 0 0.75]);
        [x_trans,y_trans]=select_multi_polygon('m');
        pset=set(htxt2,'Visible','off');
        x_trans(isnan(x_trans))=[];
        y_trans(isnan(y_trans))=[];
        n_trans=length(x_trans)/2;
        %         close(figure(12));
    end
    
end

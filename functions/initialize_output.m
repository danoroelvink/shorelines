function [O,P,V]=initialize_output(S)
% function [O,P,V]=initialize_output(S)
%
% OUTPUT
%    O       Output data structure
%    P       Output data structure projected on a grid
%    V       3D Matrix with stacked 2D frames for video
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
    
    fprintf('  Initialize output data-structures \n'); 
    O=struct;
    P=struct;
    V=struct('cdata', cell(1,1), 'colormap', cell(1,1));
       
    % prepare output on direct output
    fieldnmO = {'it','dt','tc','nt','timenum',...
                'n','x','y','distx','distQS','dSds',...
                'Wberm','xdune','ydune','qs','ql','qw','R','SWL','PHIcxy','Bf','Bm','Bfm',...
                'n1','x1','y1','distx1','distQS1','PHIc','QS',...
                'HSo','PHIo','TP','PHIf',...
                'HS','PHI','dPHI',...
                'HSbr','PHIbr','dPHIbr','hbr',...
                'x_hard','y_hard','n_hard',...
                'x_nour','y_nour','n_nour',...
                'x_fnour','y_fnour','n_fnour','V_fnour_t','q_fnour_t',... 
                'x_groyne','y_groyne'};
    for ff=1:length(fieldnmO)
        O.(fieldnmO{ff})=[];
    end
    O.outputdir=S.outputdir;
    O.storageinterval=S.storageinterval;
    O.xyprofiles=S.xyprofiles;
    % create output fields in case continuous output is needed at profiles.
    fieldnmOP={'timenum_profile','adt_profile','it_profile','c_profile','xc_profile','yc_profile','d_profile','xd_profile','yd_profile'};
    if ~isempty(S.xyprofiles)
        for ff=1:length(fieldnmOP)
            O.(fieldnmOP{ff})=[];
        end
    end
    
    % prepare output on projected grids
    if ~isempty(S.xyout)
        fieldnmP = {'it','dt','tc','nt','cntr','itout','timenum',...
                    'xg','yg','zg',...
                    'xc','yc','dist',...
                    'PHIc','PHIf',...
                    'HSo','HStdp','HSbr',...
                    'TP',...
                    'PHIo','PHItdp','PHIbr',...
                    'QS','QSmax'};
        
        if iscell(S.xyout)
            if length(S.xyout)>=1
                if length(S.xyout{1})==1 && length(S.xyout{end})==1
                    S.xyout={[cell2mat(S.xyout)]};
                end
            end
            ppval=length(S.xyout);
        else
            ppval=size(S.xyout,1);
        end
        for pp=1:ppval
            for ff=1:length(fieldnmP)
                P(pp).(fieldnmP{ff})=[];
            end

            % add reference coastline for projections
            if iscell(S.xyout)
                % use a cell with a [2x2] matrix with the startpoint of the projection grid (x1,y1) at the first row, and the end point (x2,y2) on the second row.
                % the grid will be used 1:1 and is not further interpolated
                P(pp).xg=S.xyout{pp}(:,1);
                P(pp).yg=S.xyout{pp}(:,2);
                if size(S.xyout{pp},1)==2 && size(S.xyout{pp},2)>2
                    P(pp).xg=S.xyout{pp}(1,:)';
                    P(pp).yg=S.xyout{pp}(2,:)';
                end
            else
                % use a [4x1] of (x1,y1,2,y2) matrix, with consecutively the startpoint (x1,y1) and endpoint (x2,y2) of the projection grid.
                % in this case the model will automatically make a grid in-between these points with ds0 gridsize. 
                x1=S.xyout(pp,1);
                y1=S.xyout(pp,2);
                x2=S.xyout(pp,3);
                y2=S.xyout(pp,4);
                L=hypot(x2-x1,y2-y1);
                nrcells=max(round(L./mean(S.ds0(:,end)))-1,1);
                dist1=[0,L];
                dist2=[0:L/nrcells:L];
                P(pp).xg=interp1(dist1,[x1,x2],dist2(:));
                P(pp).yg=interp1(dist1,[y1,y2],dist2(:));
            end

            % compute alongshore distance along grid
            P(pp).dist=[0;cumsum((diff(P(pp).xg).^2+diff(P(pp).yg).^2).^0.5)];

            % initialize output fields with zeros
            fields1={'TP','HSo','HStdp','HSbr','QS','QSmax'};
            fields2={'PHIc','PHIf','PHIo','PHItdp','PHIbr'};
            for ff=1:length(fields1)
            P(pp).(fields1{ff})=zeros(length(P(pp).xg),1);
            end
            for ff=1:length(fields2)
            P(pp).(fields2{ff})=zeros(length(P(pp).xg),1);
            end
            
            % reset counter for moving average of current output step of the P-structure
            P(pp).cntr=1;
            % set index for the first output step of the P-structure
            P(pp).itout=1;
        end
    end
end
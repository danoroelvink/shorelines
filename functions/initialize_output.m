function [O,P,V]=initialize_output(S,TIME,COAST,DUNE,MUD)
% function [O,P,V]=initialize_output(S,TIME,COAST,DUNE,MUD)
%
% Initialization function of the output data structures. 
%
% OUTPUT:
%    O      : Output data structure
%    P      : Output data structure projected on a grid
%    V      : 3D Matrix with stacked 2D frames for video
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
    
    fprintf('  Initialize output data-structures \n'); 
    O=struct;
    P=struct;
    V=struct('cdata', cell(1,1), 'colormap', cell(1,1));
       
    O.netcdf=S.netcdf;
    O.storageinterval=S.storageinterval;
    O.storagedate=S.storagedate;                                           % S.storagedate= {'2024-01-01', '2025-03-12'}
    O.outputdir=S.outputdir;
    if ~exist(fullfile(pwd,O.outputdir),'dir')
       mkdir(fullfile(pwd,O.outputdir));
    end
    O.outputfile=S.outputfile;
    [dirnm,filnm,extnm]=fileparts(O.outputfile);
    if O.netcdf~=1
        O.outputfile = [filnm,'.mat'];
    else
        O.outputfile = [filnm,'.nc'];
    end

    % prepare output on direct output
    O.itout=1;                                                             % counter of output timesteps
    if O.netcdf~=1
        fieldnmO = {'it','dt','tc','nt','time','timenum',...
                    'x','y','distx','h0',... % ,'dSds' %'distx1, % 'n1','x1','y1',
                    'PHIc','PHIcxy','PHIf',... 
                    'distQS','QS',...
                    'TP','HSo','HStdp','HSbr','PHIo','PHItdp','PHIbr',... %'dPHItdp','dPHIbr','hbr',... 
                    'wberm','dfelev','dcelev','xdune','ydune','qs','ql','qw','R','SWL','Bf','Bm','Bfm','hd0'... 
                    'xhard','yhard','nhard',...
                    'xnour','ynour','nnour',...
                    'x_fnour','y_fnour','n_fnour','V_fnour_t','q_fnour_t',... 
                    'x_groyne','y_groyne'};
        for ff=1:length(fieldnmO)
            O.(fieldnmO{ff})=[];
        end
    end
    
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
        fieldnmP = {'it','dt','tc','nt','cntr','itout','time','timenum',...
                   'xg','yg','distg',...
                   'zc','xc','yc','h0','volc',...
                   'PHIc','PHIcxy','PHIf','TP','HSo','HStdp','HSbr',...
                   'PHIo','PHItdp','PHIbr','QS','QSmax',...
                   'zd','xd','yd','hd0','vold','dfelev','dcelev',...
                   'wberm','qs','ql','qw','R','SWL',...
                   'Bf','Bm','Bfm'};

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
            if iscell(S.xyout) || size(S.xyout,2)==2
                % use a cell with a [2x2] matrix with the startpoint of the projection grid (x1,y1) at the first row, and the end point (x2,y2) on the second row.
                % the grid will be used 1:1 and is not further interpolated
                P(pp).xg=S.xyout{pp}(:,1);
                P(pp).yg=S.xyout{pp}(:,2);
                if size(S.xyout{pp},1)==2 && size(S.xyout{pp},2)>2
                    P(pp).xg=S.xyout{pp}(1,:)';
                    P(pp).yg=S.xyout{pp}(2,:)';
                end
            else
                % use a [4x1] of (x1,y1,x2,y2) matrix, with consecutively the startpoint (x1,y1) and endpoint (x2,y2) of the projection grid.
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
            P(pp).distg=[0;cumsum((diff(P(pp).xg).^2+diff(P(pp).yg).^2).^0.5)];

            % initialize output fields with zeros
            if O.netcdf==1
                fields={};
            else
                fields={'TP','HSo','HStdp','HSbr','QS','QSmax',...
                        'PHIc','PHIcxy','PHIf','PHIo','PHItdp','PHIbr'};
            end
            for ff=1:length(fields)
               P(pp).(fields{ff})=zeros(length(P(pp).xg),1);
            end
            
            % reset counter for moving average of current output step of the P-structure
            P(pp).cntr=1;
            % set index for the first output step of the P-structure
            P(pp).itout=1;
        end
    end

    %% NETCDF
    if O.netcdf==1
        %        USE NCname        NClongname      UNIT                    DIM  STRUC   VARfield      TYPE       
        O.Ivars={ 1, 'it'          ,'it'          ,'-'                     , 1, 'TIME'  ,'it'        ,'Instantaneous'; ...
                  1, 'dt'          ,'dt'          ,'yr'                    , 1, 'TIME'  ,'dt'        ,'Instantaneous'; ...
                  1, 'tc'          ,'tc'          ,'-'                     , 1, 'TIME'  ,'tc'        ,'Instantaneous'; ...
                  1, 'nt'          ,'nt'          ,'-'                     , 1, 'TIME'  ,'nt'        ,'Instantaneous'; ...
                  1, 'time'        ,'time'        ,'days since model start', 1, 'TIME'  ,'time'      ,'Instantaneous'; ...
                  1, 'datenum'     ,'datenum'     ,'days since Jan 01 0000', 1, 'TIME'  ,'tnow'      ,'Instantaneous'; ...
                  1, 'x'           ,'x'           ,'m'                     , 2, 'COAST' ,'x_mc'      ,'Instantaneous'; ...
                  1, 'y'           ,'y'           ,'m'                     , 2, 'COAST' ,'y_mc'      ,'Instantaneous'; ...
                  0, 'n'           ,'n'           ,'-'                     , 2, 'COAST' ,'n_mc'      ,''; ...
                  1, 'distx'       ,'distx'       ,'m'                     , 2, 'COAST' ,'distx'     ,'Instantaneous'; ...
                  1, 'distQS'      ,'distQS'      ,'m'                     , 2, 'COAST' ,'distQS'    ,'Instantaneous'; ...
                  0, 'x0'          ,'x0'          ,'m'                     , 2, 'COAST' ,'x0_mc'     ,'Instantaneous'; ...
                  0, 'y0'          ,'y0'          ,'m'                     , 2, 'COAST' ,'y0_mc'     ,'Instantaneous'; ...
                  0, 'n0'          ,'n0'          ,'-'                     , 2, 'COAST' ,'n0_mc'     ,''; ...          
                  0, 'distx0'      ,'distx0'      ,'m'                     , 2, 'COAST' ,'distx0'    ,'Instantaneous'; ...
                  0, 'distQS0'     ,'distQS0'     ,'m'                     , 2, 'COAST' ,'distQS0'   ,'Instantaneous'; ...
                  0, 'dSds'        ,'dSds'        ,'m^3/m/yr'              , 2, 'COAST' ,'dSds_mc'   ,'Instantaneous'; ...
                  1, 'QS'          ,'QS'          ,'m^3/yr'                , 2, 'COAST' ,'QS_mc'     ,'Instantaneous'; ...
                  1, 'hactive'     ,'hactive'     ,'m'                     , 2, 'COAST' ,'h0_mc'     ,'Instantaneous'; ...
                  1, 'PHIc'        ,'PHIc'        ,'°N'                    , 2, 'COAST' ,'PHIc'      ,''; ...
                  1, 'PHIcxy'      ,'PHIcxy'      ,'°N'                    , 2, 'COAST' ,'PHIcxy'    ,''; ...
                  1, 'PHIf'        ,'PHIf'        ,'°N'                    , 2, 'COAST' ,'PHIf_mc'   ,'Instantaneous'; ...
                  1, 'TP'          ,'TP'          ,'s'                     , 2, 'WAVE'  ,'TP_mc'     ,'Instantaneous'; ...
                  1, 'HSo'         ,'HSo'         ,'m'                     , 2, 'WAVE'  ,'HSo_mc'    ,'Instantaneous'; ...
                  1, 'HStdp'       ,'HStdp'       ,'m'                     , 2, 'WAVE'  ,'HStdp_mc'  ,'Instantaneous'; ...
                  1, 'HSbr'        ,'HSbr'        ,'m'                     , 2, 'WAVE'  ,'HSbr_mc'   ,'Instantaneous'; ...
                  1, 'PHIo'        ,'PHIo'        ,'°N'                    , 2, 'WAVE'  ,'PHIo_mc'   ,'Instantaneous'; ...
                  1, 'PHItdp'      ,'PHItdp'      ,'°N'                    , 2, 'WAVE'  ,'PHItdp_mc' ,'Instantaneous'; ...
                  1, 'PHIbr'       ,'PHIbr'       ,'°N'                    , 2, 'WAVE'  ,'PHIbr_mc'  ,'Instantaneous'; ...
                  0, 'dPHIo'       ,'dPHIo'       ,'°'                     , 2, 'WAVE'  ,'dPHIo_mc'  ,''; ...
                  0, 'dPHItdp'     ,'dPHItdp'     ,'°'                     , 2, 'WAVE'  ,'dPHItdp_mc','Instantaneous'; ...
                  0, 'dPHIbr'      ,'dPHIbr'      ,'°'                     , 2, 'WAVE'  ,'dPHIbr_mc' ,'Instantaneous'; ...
                  0, 'hbr'         ,'hbr'         ,'m'                     , 2, 'WAVE'  ,'hbr_mc'    ,'Instantaneous'; ...
                  0, 'xhard'       ,'xhard'       ,'m'                     , 1, 'STRUC' ,'xhard'     ,''; ...
                  0, 'yhard'       ,'yhard'       ,'m'                     , 1, 'STRUC' ,'yhard'     ,''; ...
                  0, 'nhard'       ,'nhard'       ,'-'                     , 1, 'STRUC' ,'nhard'     ,''; ...
                  0, 'xrevet'      ,'xrevet'      ,'m'                     , 1, 'STRUC' ,'xrevet'    ,''; ...
                  0, 'yrevet'      ,'yrevet'      ,'m'                     , 1, 'STRUC' ,'yrevet'    ,''; ...
                  0, 'xnour'       ,'xnour'       ,'m'                     , 1, 'STRUC' ,'xnour'     ,''; ...
                  0, 'ynour'       ,'ynour'       ,'m'                     , 1, 'STRUC' ,'ynour'     ,''; ...
                  0, 'nnour'       ,'nnour'       ,'-'                     , 2, 'STRUC' ,'nnour'     ,''; ...
                  0, 'tnour'       ,'tnour'       ,'days since Jan 01 0000', 2, 'STRUC' ,'tnour'     ,''; ...
                  0, 'vnour'       ,'vnour'       ,'m^3'                   , 2, 'STRUC' ,'vnour'     ,''; ...
                  0, 'x_fnour'     ,'x_fnour'     ,'m'                     , 1, 'STRUC' ,'x_fnour'   ,''; ...
                  0, 'y_fnour'     ,'y_fnour'     ,'m'                     , 1, 'STRUC' ,'y_fnour'   ,''; ...
                  0, 'n_fnour'     ,'n_fnour'     ,'m'                     , 1, 'STRUC' ,'n_fnour'   ,''; ...
                  0, 't_fnour'     ,'t_fnour'     ,'days since Jan 01 0000', 1, 'STRUC' ,'t_fnour'   ,''; ...
                  0, 'V_fnour_t'   ,'V_fnour_t'   ,'m^3'                   , 1, 'STRUC' ,'V_fnour_t' ,''; ... 
                  0, 'q_fnour_t'   ,'q_fnour_t'   ,'m^3/m'                 , 1, 'STRUC' ,'q_fnour_t' ,''; ...         
                  1, 'x_dune'      ,'x_dune'      ,'m'                     , 2, 'COAST' ,'xdune_mc'  ,''; ...
                  1, 'y_dune'      ,'y_dune'      ,'m'                     , 2, 'COAST' ,'ydune_mc'  ,''; ...
                  1, 'wberm_dune'  ,'wberm_dune'  ,'m'                     , 2, 'COAST' ,'wberm_mc'  ,''; ...
                  1, 'hactive_dune','hactive_dune','m'                     , 2, 'COAST' ,'h0dune_mc' ,''; ...
                  1, 'dfelev_dune' ,'dfelev_dune' ,'m'                     , 2, 'COAST' ,'dfelev_mc' ,''; ...
                  1, 'dcelev_dune' ,'dcelev_dune' ,'m'                     , 2, 'COAST' ,'dcelev_mc' ,''; ...
                  1, 'qs_dune'     ,'qs_dune'     ,'m^3/m/yr'              , 2, 'COAST' ,'qs_mc'     ,''; ...
                  1, 'ql_dune'     ,'ql_dune'     ,'m^3/m/yr'              , 2, 'COAST' ,'ql_mc'     ,''; ...
                  1, 'qw_dune'     ,'qw_dune'     ,'m^3/m/yr'              , 2, 'COAST' ,'qw_mc'     ,''; ...
                  1, 'R_dune'      ,'R_dune'      ,'m'                     , 2, 'COAST' ,'R_mc'      ,''; ...
                  1, 'SWL_dune'    ,'SWL_dune'    ,'m'                     , 2, 'COAST' ,'SWL_mc'    ,''; ...
                  1, 'Bf_mud'      ,'Bf_mud'      ,'m'                     , 2, 'COAST' ,'Bf'        ,''; ...
                  1, 'Bm_mud'      ,'Bm_mud'      ,'m'                     , 2, 'COAST' ,'Bm'        ,''; ...
                  1, 'Bfm_mud'     ,'Bfm_mud'     ,'m'                     , 2, 'COAST' ,'Bfm'       ,''  };
               
        O.Pvars={ 1, 'it'          ,'it'          ,'-'                     , 1,'TIME'  ,'it'        , 'Instantaneous'; ...
                  1, 'dt'          ,'dt'          ,'yr'                    , 1,'TIME'  ,'dt'        , 'Instantaneous'; ...
                  1, 'tc'          ,'tc'          ,'-'                     , 1,'TIME'  ,'tc'        , 'Instantaneous'; ...
                  1, 'nt'          ,'nt'          ,'-'                     , 1,'TIME'  ,'nt'        , 'Instantaneous'; ...
                  1, 'time'        ,'time'        ,'days since model start', 1,'TIME'  ,'time'      , 'Instantaneous'; ...
                  1, 'datenum'     ,'datenum'     ,'days since Jan 01 0000', 1,'TIME'  ,'timenum'   , 'Instantaneous'; ...
                  1, 'xg'          ,'xg'          ,'m'                     , 1,'COAST' ,'xg'        , 'Projected'; ...
                  1, 'yg'          ,'yg'          ,'m'                     , 1,'COAST' ,'yg'        , 'Projected'; ...          
                  1, 'distg'       ,'distg'       ,'m'                     , 2,'COAST' ,'distg'     , 'Projected'; ...          
                  1, 'zc'          ,'zc'          ,'m'                     , 2,'COAST' ,'zc'        , 'Projected'; ...
                  1, 'xc'          ,'xc'          ,'m'                     , 2,'COAST' ,'xc'        , 'Projected'; ...
                  1, 'yc'          ,'yc'          ,'m'                     , 2,'COAST' ,'yc'        , 'Projected'; ...          
                  1, 'hactive'     ,'hactive'     ,'m'                     , 2,'COAST' ,'h0'        , 'Projected'; ...          
                  1, 'volc'        ,'volc'        ,'m^3/m'                 , 2,'COAST' ,'volc'      , 'Projected'; ...          
                  1, 'PHIc'        ,'PHIc'        ,'°N'                    , 2,'COAST' ,'PHIc'      , 'Projected'; ...
                  1, 'PHIcxy'      ,'PHIcxy'      ,'°N'                    , 2,'COAST' ,'PHIcxy'    , 'Projected'; ...
                  1, 'PHIf'        ,'PHIf'        ,'°N'                    , 2,'COAST' ,'PHIf'      , 'Projected'; ...
                  1, 'TP'          ,'TP'          ,'s'                     , 2,'COAST' ,'TP'        , 'Projected'; ...
                  1, 'HSo'         ,'HSo'         ,'m'                     , 2,'COAST' ,'HSo'       , 'Projected'; ...
                  1, 'HStdp'       ,'HStdp'       ,'m'                     , 2,'COAST' ,'HStdp'     , 'Projected'; ...
                  1, 'HSbr'        ,'HSbr'        ,'m'                     , 2,'COAST' ,'HSbr'      , 'Projected'; ...
                  1, 'PHIo'        ,'PHIo'        ,'°N'                    , 2,'COAST' ,'PHIo'      , 'Projected'; ...
                  1, 'PHItdp'      ,'PHItdp'      ,'°N'                    , 2,'COAST' ,'PHItdp'    , 'Projected'; ...
                  1, 'PHIbr'       ,'PHIbr'       ,'°N'                    , 2,'COAST' ,'PHIbr'     , 'Projected'; ...
                  1, 'QS'          ,'QS'          ,'m^3/yr'                , 2,'COAST' ,'QS'        , 'Projected'; ...
                  1, 'QSmax'       ,'QSmax'       ,'m^3/yr'                , 2,'COAST' ,'QSmax'     , 'Projected'; ...
                  1, 'zd_dune'     ,'zd_dune'     ,'m'                     , 2,'COAST' ,'zd'        , 'Projected'; ...
                  1, 'xd_dune'     ,'xd_dune'     ,'m'                     , 2,'COAST' ,'xd'        , 'Projected'; ...
                  1, 'yd_dune'     ,'yd_dune'     ,'m'                     , 2,'COAST' ,'yd'        , 'Projected'; ...
                  1, 'hactive_dune','hactive_dune','m'                     , 2,'COAST' ,'hd0'       , 'Projected'; ...
                  1, 'vol_dune'    ,'vol_dune'    ,'m^3/m'                 , 2,'COAST' ,'vold'      , 'Projected'; ...
                  1, 'dfelev_dune' ,'dfelev_dune' ,'m'                     , 2,'COAST' ,'dfelev'    , 'Projected'; ...
                  1, 'dcelev_dune' ,'dcelev_dune' ,'m'                     , 2,'COAST' ,'dcelev'    , 'Projected'; ...
                  1, 'wberm_dune'  ,'wberm_dune'  ,'m'                     , 2,'COAST' ,'wberm'     , 'Projected'; ...
                  1, 'qs_dune'     ,'qs_dune'     ,'m^3/m/yr'              , 2,'COAST' ,'qs'        , 'Projected'; ...
                  1, 'ql_dune'     ,'ql_dune'     ,'m^3/m/yr'              , 2,'COAST' ,'ql'        , 'Projected'; ...
                  1, 'qw_dune'     ,'qw_dune'     ,'m^3/m/yr'              , 2,'COAST' ,'qw'        , 'Projected'; ...
                  1, 'R_dune'      ,'R_dune'      ,'m'                     , 2,'COAST' ,'R'         , 'Projected'; ...
                  1, 'SWL_dune'    ,'SWL_dune'    ,'m'                     , 2,'COAST' ,'SWL'       , 'Projected'; ...
                  1, 'Bf_mud'      ,'Bf_mud'      ,'m'                     , 2,'COAST' ,'Bf'        , 'Projected'; ...
                  1, 'Bm_mud'      ,'Bm_mud'      ,'m'                     , 2,'COAST' ,'Bm'        , 'Projected'; ...
                  1, 'Bfm_mud'     ,'Bfm_mud'     ,'m'                     , 2,'COAST' ,'Bfm'       , 'Projected'  };

        % Optimized NetCDF writer with unlimited 'space' dimension
        O.nPspace = isempty(fieldnames(P)) * 0 + ~isempty(fieldnames(P)) * numel(P);
        O.ncinitialized=[];
        O.separatepgrids=S.separatepgrids;

        % Check init of ncfile
        if isempty(O.ncinitialized)
            %% === 1) CREATE & DEFINE FILE ===
            ns = numel(COAST.x_mc);
            O.ncfile = fullfile(pwd,O.outputdir,O.outputfile);
            [P]=initialize_netcdf(O, P, ns, DUNE.used, MUD.used);
            
            %% === 2) OPEN ONCE & DISABLE FILL ===
            O.ncid = netcdf.open(O.ncfile, 'NC_WRITE');
            
            netcdf.reDef(O.ncid);
            netcdf.setFill(O.ncid,'NC_FILL');
            netcdf.endDef(O.ncid);
            
            %% === 3) CACHE VARIABLE IDs ===
            % time?scalars
            O.ncIDs = zeros(1,size(O.Ivars,1));
            for i = 1:size(O.Ivars,1)
                % switch for using/exporting this variable + check if for DUNE or MUD variables
                if O.Ivars{i,1}==1 && (isempty(findstr(O.Ivars{i,2},'dune')) || DUNE.used==1) && (isempty(findstr(O.Ivars{i,2},'mud')) || MUD.used==1)  
                    O.ncIDs(i) = netcdf.inqVarID(O.ncid, O.Ivars{i,2});
                end
            end
            
            % projected structures P
            O.ncPInstIDs = zeros(O.nPspace,size(O.Pvars,1));
            for pp = 1:O.nPspace
                
                if O.separatepgrids==1
                    % open nc files of p-grids
                    P(pp).ncid = netcdf.open(P(pp).ncfile, 'NC_WRITE');
                    netcdf.reDef(P(pp).ncid);
                    netcdf.setFill(P(pp).ncid,'NC_FILL');
                    netcdf.endDef(P(pp).ncid);
                else
                    % get group
                    grp = sprintf('projected_grid_%03d',pp);
                    P(pp).ncid = netcdf.inqNcid(O.ncid, grp);
                end
                
                for i = 1:size(O.Pvars,1)
                   % switch for using/exporting this variable + check if for DUNE or MUD variables
                   if O.Pvars{i,1}==1 && (isempty(findstr(O.Pvars{i,2},'dune')) || DUNE.used==1) && (isempty(findstr(O.Pvars{i,2},'mud')) || MUD.used==1)  
                      O.ncPInstIDs(pp,i) = netcdf.inqVarID(P(pp).ncid, O.Pvars{i,2});
                   end
                end
                if O.separatepgrids==1
                   netcdf.close(P(pp).ncid);
                   P(pp).ncinitialized = true;
                end
            end
            
            netcdf.close(O.ncid);
            O.ncinitialized = true;

        end
    end 
end

%% ---------------------------------------------------------------------
function [P]=initialize_netcdf(O, P, ns, save_dune, save_mud)

   if exist(O.ncfile,'file') == 2
      delete(O.ncfile)
   end

   deflateLevel = 5;
   cmode = bitor( netcdf.getConstant('CLOBBER'), ...
                  netcdf.getConstant('NETCDF4') );
   ncid = netcdf.create(O.ncfile, cmode);

   %--- dimensions ---
   time_dim  = netcdf.defDim(ncid,'time',  netcdf.getConstant('NC_UNLIMITED'));
   space_dim = netcdf.defDim(ncid,'space', netcdf.getConstant('NC_UNLIMITED'));
   grid_dim  = netcdf.defDim(ncid,'grid',  netcdf.getConstant('NC_UNLIMITED'));

   %--- time scalar variables ---     
   for v = 1:size(O.Ivars,1)
      % Example line from the table
      %  USE NCname        NClongname      UNIT                    DIM  STRUC   VARfield      TYPE       
      %   1, 'it'          ,'it'          ,'-'                     , 1, 'TIME'  ,'it'        ,'Instantaneous'; ...
      
      % switch for using/exporting this variable + check if for DUNE or MUD variables
      if O.Ivars{v,1}==1 && (isempty(strfind(O.Ivars{v,2},'dune')) || save_dune==1) && (isempty(strfind(O.Ivars{v,2},'mud')) || save_mud==1)        

         if O.Ivars{v,5}==1   % switch for 1D or 2D variable (1/2)
            varid = netcdf.defVar(ncid, O.Ivars{v,2}, 'NC_DOUBLE', time_dim);
            netcdf.defVarDeflate(ncid, varid, true, true, deflateLevel);
            netcdf.defVarChunking(ncid, varid, 'CHUNKED', 1);
         else
            varid = netcdf.defVar(ncid, O.Ivars{v,2}, 'NC_DOUBLE', [space_dim, time_dim]);
            netcdf.defVarDeflate(ncid, varid, true, true, deflateLevel);
            netcdf.defVarChunking(ncid, varid, 'CHUNKED', [min(1000,ns),1]);
         end
         netcdf.putAtt(ncid, varid, 'long_name', O.Ivars{v,3});
         netcdf.putAtt(ncid, varid, 'units',     O.Ivars{v,4});
         netcdf.defVarFill(ncid, varid, false, NaN);
      end
   end
   
   %--- projected P-structures (fixed grid time) ---
   if ~isempty(fields(P))
      for pp = 1:length(P) 
          
         if O.separatepgrids==1
            % store in seprate ncfiles
            P(pp).ncfile = [O.ncfile(1:end-3),'_grid',num2str(pp,'%02.0f'),'.nc'];
            if exist(P(pp).ncfile,'file') == 2
               delete(P(pp).ncfile)
            end
            ncidgrp = netcdf.create(P(pp).ncfile, cmode);
    
            %--- dimensions ---
            time_dim  = netcdf.defDim(ncidgrp,'time',  netcdf.getConstant('NC_UNLIMITED'));
            space_dim = netcdf.defDim(ncidgrp,'space', netcdf.getConstant('NC_UNLIMITED'));
            grid_dim  = netcdf.defDim(ncidgrp,'grid',  netcdf.getConstant('NC_UNLIMITED'));
         else
            % open as sub-group in ncfile
            ncidgrp = netcdf.defGrp(ncid,sprintf('projected_grid_%03d',pp));
            np      = numel(P(pp).xg);
            time_dim = netcdf.defDim(ncidgrp,'time',  netcdf.getConstant('NC_UNLIMITED'));
            space_dim= netcdf.defDim(ncidgrp,'space', np);
         end

         for v = 1:size(O.Pvars,1)
            % switch for using/exporting this variable + check if for DUNE or MUD variables
            if O.Pvars{v,1}==1 && (isempty(strfind(O.Pvars{v,2},'dune')) || save_dune==1) && (isempty(strfind(O.Pvars{v,2},'mud')) || save_mud==1)  
               if strcmpi(O.Pvars{v,2},'xg') || strcmpi(O.Pvars{v,2},'yg')  % switch for 1D or 2D variable (1/2)
                  varid = netcdf.defVar(ncidgrp, O.Pvars{v,2}, 'NC_DOUBLE', space_dim);
                  netcdf.defVarDeflate(ncidgrp, varid, true, true, deflateLevel);
                  netcdf.defVarChunking(ncidgrp, varid, 'CHUNKED', space_dim);
               elseif O.Pvars{v,5}==1   % switch for 1D or 2D variable (1/2)
                  varid = netcdf.defVar(ncidgrp, O.Pvars{v,2}, 'NC_DOUBLE', time_dim);
                  netcdf.defVarDeflate(ncidgrp, varid, true, true, deflateLevel);
                  netcdf.defVarChunking(ncidgrp, varid, 'CHUNKED', 1);
               else
                  varid = netcdf.defVar(ncidgrp, O.Pvars{v,2}, 'NC_DOUBLE', [space_dim, time_dim]);
                  netcdf.defVarDeflate(ncidgrp, varid, true, true, deflateLevel);
                  netcdf.defVarChunking(ncidgrp, varid, 'CHUNKED', [min(1000,space_dim),1]);
               end
               netcdf.putAtt(ncidgrp, varid, 'long_name', O.Pvars{v,3});
               netcdf.putAtt(ncidgrp, varid, 'units',     O.Pvars{v,4});
               netcdf.defVarFill(ncidgrp, varid, false, NaN);
            end
         end

         if O.separatepgrids==1
           %--- global metadata ---
            netcdf.putAtt(ncidgrp, netcdf.getConstant('NC_GLOBAL'), 'title', 'ShorelineS');
            netcdf.putAtt(ncidgrp, netcdf.getConstant('NC_GLOBAL'), 'institution','IHE Delft/Deltares');
            netcdf.putAtt(ncidgrp, netcdf.getConstant('NC_GLOBAL'), 'history', ['Created on ' datestr(now,'yyyy-mm-dd HH:MM:SS')]);
            netcdf.putAtt(ncidgrp, netcdf.getConstant('NC_GLOBAL'), 'references', 'Roelvink et al. (2020). https://doi.org/10.3389/fmars.2020.00535');
            netcdf.putAtt(ncidgrp, netcdf.getConstant('NC_GLOBAL'), 'Conventions', 'CF-1.8');
        
            netcdf.endDef(ncidgrp);
            netcdf.close(ncidgrp);
         else
            netcdf.endDef(ncidgrp);
         end
      end
   end

   %--- global metadata ---
   netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', 'ShorelineS');
   netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'institution','IHE Delft/Deltares');
   netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'history', ['Created on ' datestr(now,'yyyy-mm-dd HH:MM:SS')]);
   netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'references', 'Roelvink et al. (2020). https://doi.org/10.3389/fmars.2020.00535');
   netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'Conventions', 'CF-1.8');

   netcdf.endDef(ncid);
   netcdf.close(ncid);
end
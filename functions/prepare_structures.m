function [STRUC]=prepare_structures(S,COAST)
% function [STRUC]=prepare_structures(S,COAST)
%
% Initializes the structures used in the model (e.g. revetments, groynes and offshore breakwaters).
% Also, the parameters for wave transmission at offshore breakwaters are read and prepared. 
% The relevant input information is stored in the data-structure STRUC. 
%
% INPUT:
%     S
%         .ldbstructures     : name of file with groynes and offshore breakwaters (ASCII, [Nx2] without header)
%         .struc             : switch for using structures (0/1)
%         .xhard             : x-coordinates of groynes and offshore breakwaters [m]
%         .yhard             : y-coordinates of groynes and offshore breakwaters [m]
%         .diffraction       : switch for using diffraction (0/1)
%         .diffdist          : critical distance for sheltering of diffracted waves [m] (i.e. relevant for bays)
%         .kdform            : diffraction approach for the wave height reduction (either 'Kamphuis' and 'Roelvink')
%         .wdform            : diffraction approach for the directional spreading of the waves (either 'Roelvink' of 'Dabees')
%         .transmission      : switch for using transmission (0/1)
%         .transmform        : used formulation for calculating transmission 
%         .transmbwdepth     : depth at breakwater location [m]
%         .transmcrestheight : breakwater crest height [m]
%         .transmslope       : breakwater slope [-]
%         .transmcrestwidth  : breakwater crest width [m]
%         .transmdir         : switch for using either the incoming wave direction (0) or a weighted direction of the incoming wave and shore-normal of the structure (1)
%         .transmd50         : D50 of breakwater armour material [m]
%         .ldbrevetments     : name of file with revetments (ASCII, [Nx2] without header)
%         .revet             : switch for using revetments (0/1)
%         .xrevet            : x-coordinates of revetments [m]
%         .yrevet            : y-coordinates of revetments [m]
%         .iterrev           : iterations used to determine alongshore transport along the revetment
%         .ldbpermeable      : name of file with permeable structures (ASCII, [Nx2] without header)
%         .perm              : switch for using permeable structures (0/1)
%         .xperm             : x-coordinates of permeable structures [m]
%         .yperm             : y-coordinates of permeable structures [m]
%         .xyoffset          : offset of the x and y coordinates (to make numbers on axis smaller) [1x2]
%
% OUTPUT:
%     STRUC
%         .xhard             : x-coordinates of groynes and offshore breakwaters [m]
%         .yhard             : y-coordinates of groynes and offshore breakwaters [m]
%         .nhard             : number of groynes and offshore breakwaters
%         .xrevet            : x-coordinates of revetments [m]
%         .yrevet            : y-coordinates of revetments [m]
%         .nrevet            : number of revetments
%         .xperm             : x-coordinates of permeable structures [m]
%         .yperm             : y-coordinates of permeable structures [m]
%         .nperm             : number of permeable structures
%         .transmission      : switch for using transmission (0/1)
%         .transmform	     : used formulation for calculating transmission 
%         .transmbwdepth     : depth at breakwater location [m]
%         .transmcrestheight : breakwater crest height [m]
%         .transmslope	     : breakwater slope [-]
%         .transmcrestwidth  : breakwater crest width [m]
%         .transmdir	     : switch for using either the incoming wave direction (0) or a weighted direction of the incoming wave and shore-normal of the structure (1)
%         .transmd50	     : D50 of breakwater armour material [m]
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

    fprintf('  Prepare structures \n');
    STRUC=struct;
    
    % Wave blocking structure
    STRUC.xhard=[];
    STRUC.yhard=[];
    STRUC.shard=[];
    STRUC.nhard=0;
    STRUC.type=S.structtype;
    STRUC.transmission=S.transmission;
    STRUC.transmform=S.transmform;
    STRUC.transmdir=S.transmdir;
    STRUC.bypasscontrfac=S.bypasscontrfac;                                            % the maximum transport when the coastline is at the end of the structure (always >=1). Setting the bypass fraction larger than 1 means that the accretion does not go to the tip of the structure. 
    
    % Wave transmitting structures
    STRUC.perm=S.perm;                                                                  % switch for using hard structures
    STRUC.ldbpermeable=S.ldbpermeable;                                                              % LDB with perm structures ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
    STRUC.xperm=[];     % for permeable structures
    STRUC.yperm=[];     % for permeable structures
    STRUC.wavetransm=S.wavetransm;
    
    % Revetments     
    STRUC.revet=S.revet;
    STRUC.xrevet=[];     % for permeable structures
    STRUC.yrevet=[];     % for permeable structures
    STRUC.iterrev=S.iterrev; % iterations used to determine alongshore transport along the revetment
    
    % Diffraction at a structure
    STRUC.diffraction=S.diffraction;
    STRUC.rotfac=S.rotfac;
    STRUC.kdform=S.kdform;
    STRUC.wdform=S.wdform;
    STRUC.xtip=[];       % for diffraction
    STRUC.ytip=[];       % for diffraction
    STRUC.hstip=[];      % for diffraction
    STRUC.delta0=[];     % for diffraction
    STRUC.tip_wet=[];    % for diffraction
    STRUC.xp=[];         % for diffraction
    STRUC.yp=[];         % for diffraction
    STRUC.wetstr_mc=[];  % for diffraction
    STRUC.diffdist=S.diffdist;
    
    % Groyne
    STRUC.bypassdistpwr=S.bypassdistpwr;
    STRUC.groinelev=S.groinelev;
    
    %% Coastal structures (blocking waves)
    if S.struct || ~isempty(S.ldbstructures)
        %  add groynes/breakwaters interactively with the graphical user interface
        if strcmpi(S.ldbstructures,'manual') || strcmpi(S.ldbstructures,'interactive') 
            figure(11);
            axis equal;
            xl=xlim;yl=ylim;
            htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add structure (LMB); Next structure (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
            [xhard,yhard]=select_multi_polygon('k');
            set(htxt2,'Visible','off');
            STRUC.xhard=xhard;
            STRUC.yhard=yhard;

        %  try to load a file with groyne/breakwater 
        elseif ~isempty(S.ldbstructures)
            xyhard=load(S.ldbstructures);
            STRUC.xhard=xyhard(:,1)'-S.xyoffset(1);
            STRUC.yhard=xyhard(:,2)'-S.xyoffset(2);

        % add a groyne by directly specifying the xhard and yhard
        elseif ~isempty(S.xhard) && ~isempty(S.yhard)
            STRUC.xhard=S.xhard(:)'-S.xyoffset(1);
            STRUC.yhard=S.yhard(:)'-S.xyoffset(2);
        % no groynes or breakwaters
        else
            S.struct=0;
            STRUC.xhard=[];
            STRUC.yhard=[];
        end
    end
    if ~isempty(STRUC.xhard)
        STRUC.nhard=length(find(isnan(STRUC.xhard)))+1;
    end
    
    if length(S.bypasscontrfac)<STRUC.nhard
        STRUC.bypasscontrfac(1:STRUC.nhard)=S.bypasscontrfac(1);
    elseif length(S.bypasscontrfac)==STRUC.nhard
        STRUC.bypasscontrfac=S.bypasscontrfac;
    end
    
    if STRUC.transmission == 1
        if ~isempty(S.transmfile)
            charac=load(S.transmfile);
            if size(charac(1,:),2)<STRUC.nhard
                STRUC.transmbwdepth(1:STRUC.nhard)=charac(1,:);
                STRUC.transmcrestheight(1:STRUC.nhard)=charac(2,:);
                STRUC.transmslope(1:STRUC.nhard)=charac(3,:);
                STRUC.transmcrestwidth(1:STRUC.nhard)=charac(4,:);
                if length(charac(:,1)) == 4
                    STRUC.transmd50(1:STRUC.nhard)=S.transmd50; % default value
                else
                    STRUC.transmd50(1:STRUC.nhard)=charac(5,:);
                end
            elseif size(charac(1,:),2)==STRUC.nhard
                STRUC.transmbwdepth=charac(1,:);                           
                STRUC.transmcrestheight=charac(2,:);
                STRUC.transmslope=charac(3,:);
                STRUC.transmcrestwidth=charac(4,:);
                STRUC.transmd50=charac(5,:);
            end
        elseif ~isempty(STRUC.xhard) && ~isempty(STRUC.yhard)
            if length(S.transmbwdepth)<STRUC.nhard
                STRUC.transmbwdepth(1:STRUC.nhard)=S.transmbwdepth(1);
            elseif length(S.transmbwdepth)==STRUC.nhard
               STRUC.transmbwdepth=S.transmbwdepth;
            end
            if length(S.transmcrestheight)<STRUC.nhard
                STRUC.transmcrestheight(1:STRUC.nhard)=S.transmcrestheight(1);
            elseif length(S.transmcrestheight)==STRUC.nhard
               STRUC.transmcrestheight=S.transmcrestheight;
            end
            if length(S.transmslope)<STRUC.nhard
               STRUC.transmslope(1:STRUC.nhard)=S.transmslope(1);
            elseif length(S.transmslope)==STRUC.nhard
               STRUC.transmslope=S.transmslope;
            end
            if length(S.transmcrestwidth)<STRUC.nhard
               STRUC.transmcrestwidth(1:STRUC.nhard)=S.transmcrestwidth(1);
            elseif length(S.transmcrestwidth)==STRUC.nhard
               STRUC.transmcrestwidth=S.transmcrestwidth;
            end
            if length(S.transmd50)<STRUC.nhard
                STRUC.transmd50(1:STRUC.nhard)=S.transmd50(1);
            elseif length(S.transmd50)==STRUC.nhard
                STRUC.transmd50=S.transmd50;
            end 
        end
        if isfield(S,'transmform')
            STRUC.transmform=S.transmform;
        else
            STRUC.transmform='angr';
        end        
        if ischar(STRUC.transmform)
            STRUC.transmform={STRUC.transmform};
        end
        if length(STRUC.transmform)==1
            STRUC.transmform=repmat(STRUC.transmform,[1,STRUC.nhard]);
        end
    else 
        STRUC.transmbwdepth=[];
        STRUC.transmcrestheight=[];
        STRUC.transmslope=[];
        STRUC.transmcrestwidth=[];
        STRUC.transmd50=[];
    end 

    %% Permeable structures
    if S.perm 
        %  add permeable structures interactively with the graphical user interface
        if strcmpi(S.ldbpermeable,'manual') || strcmpi(S.ldbpermeable,'interactive') 
            figure(11);
            axis equal;
            xl=xlim;yl=ylim;
            htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add permeable structure (LMB); Next structure (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
            [xperm,yperm]=select_multi_polygon('k');
            set(htxt2,'Visible','off');
            STRUC.xperm=xperm;
            STRUC.yperm=yperm;
        
        %  try to load a file with permeable structures 
        elseif ~isempty(S.ldbpermeable)
            xyperm=load(S.ldbpermeable);
            STRUC.xperm=xyperm(:,1)'-S.xyoffset(1);
            STRUC.yperm=xyperm(:,2)'-S.xyoffset(2);

        % add a permeable structure by directly specifying the xperm and yperm
        elseif ~isempty(S.xperm) && ~isempty(S.yperm)
            STRUC.xperm=S.xperm(:)'-S.xyoffset(1);
            STRUC.yperm=S.yperm(:)'-S.xyoffset(2);

        % no permeable structures
        else
            S.perm=0;
            STRUC.xperm=[];
            STRUC.yperm=[];
        end
        STRUC.nperm=min(sum(isnan(STRUC.xperm))+1,length(STRUC.xperm));
    end
    
    %% Revetments
    STRUC.nrevet=0;
    if S.revet 
        %  add revetments interactively with the graphical user interface
        if strcmpi(S.ldbrevetments,'manual') || strcmpi(S.ldbrevetments,'interactive') 
            figure(11);
            axis equal;
            xl=xlim;yl=ylim;
            htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add revetment (LMB); Next revetment (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
            [xrevet,yrevet]=select_multi_polygon('k');
            set(htxt2,'Visible','off');
            STRUC.xrevet=xrevet;
            STRUC.yrevet=yrevet;

        %  try to load a file with revetment 
        elseif ~isempty(S.ldbrevetments)
            xyrevet=load(S.ldbrevetments);
            STRUC.xrevet=xyrevet(:,1)'-S.xyoffset(1);
            STRUC.yrevet=xyrevet(:,2)'-S.xyoffset(2);

        % add a revetment by directly specifying the xrevet and yrevet
        elseif ~isempty(S.xrevet) && ~isempty(S.yrevet)
            STRUC.xrevet=S.xrevet(:)'-S.xyoffset(1);
            STRUC.yrevet=S.yrevet(:)'-S.xyoffset(2);

        % no revetments
        else
            S.revet=0;
            STRUC.xrevet=[];
            STRUC.yrevet=[];
        end
    end
    if ~isempty(STRUC.xrevet) 
        % STRUC.nrevet=sum(isnan(STRUC.xrevet))+1;
        % STRUC.xhard=[STRUC.xhard,nan,STRUC.xrevet];
        % STRUC.yhard=[STRUC.yhard,nan,STRUC.yrevet];
    end
    
    % Check whether to include diffraction
    STRUC.diffraction=S.diffraction;
    if isempty(STRUC.xhard) && isempty(STRUC.xrevet) && isempty(S.xsedlim) && isempty(S.ldbsedlim) 
        STRUC.diffraction=0; 
    end
    
    % tidy up the structures without nans at the start and end
    idnotnan=find(~isnan(STRUC.xhard));
    if ~isempty(STRUC.xhard)
        STRUC.xhard=STRUC.xhard(idnotnan(1):idnotnan(end));
        STRUC.yhard=STRUC.yhard(idnotnan(1):idnotnan(end));
        idnan=find(isnan(STRUC.xhard));
        iduse=setdiff([1:length(STRUC.xhard)],idnan(diff(idnan)==1));
        STRUC.xhard=STRUC.xhard(iduse);
        STRUC.yhard=STRUC.yhard(iduse);
    end
    
    % tidy up the permeable structures without nans at the start and end
    idnotnan=find(~isnan(STRUC.xperm));
    if ~isempty(STRUC.xperm)
        STRUC.xperm=STRUC.xperm(idnotnan(1):idnotnan(end));
        STRUC.yperm=STRUC.yperm(idnotnan(1):idnotnan(end));
        idnan=find(isnan(STRUC.xperm));
        iduse=setdiff([1:length(STRUC.xperm)],idnan(diff(idnan)==1));
        STRUC.xperm=STRUC.xperm(iduse);
        STRUC.yperm=STRUC.yperm(iduse);
    end

    % tidy up the revetments without nans at the start and end
    idnotnan=find(~isnan(STRUC.xrevet));
    if ~isempty(STRUC.xrevet)
        STRUC.xrevet=STRUC.xrevet(idnotnan(1):idnotnan(end));
        STRUC.yrevet=STRUC.yrevet(idnotnan(1):idnotnan(end));
        idnan=find(isnan(STRUC.xrevet));
        iduse=setdiff([1:length(STRUC.xrevet)],idnan(diff(idnan)==1));
        STRUC.xrevet=STRUC.xrevet(iduse);
        STRUC.yrevet=STRUC.yrevet(iduse);
    end
       
    %% IDENTIFY THE STRUCTURES WHICH ARE LOCATED IN THE SEA (OR AT LAND)
    nmc=length(COAST.x_mc)-1;
    for ist=1:length(STRUC.xhard)
        [~,icl]=min(hypot(COAST.x_mc-STRUC.xhard(ist),COAST.y_mc-STRUC.yhard(ist)));
        if ~isnan(STRUC.xhard(ist))
            im1=max(icl-1,1);
            ip1=min(icl+1,nmc+1);
            if isnan(COAST.x_mc(im1)); im1=im1-1;if im1==0;im1=2;end; end
            if isnan(COAST.x_mc(ip1)); ip1=ip1+1;if ip1>nmc+1;ip1=nmc;end; end
            dirm=360-atan2d(COAST.y_mc(ip1)-COAST.y_mc(im1),COAST.x_mc(ip1)-COAST.x_mc(im1)); 
            dirstr=atan2d(STRUC.xhard(ist)-COAST.x_mc(icl),STRUC.yhard(ist)-COAST.y_mc(icl));
            STRUC.wetstr_mc(ist)=int8(cosd(dirstr-dirm)>0);  % for octave, logicals cannot be assigned NaNs later on
        else
            STRUC.wetstr_mc(ist)=int8(0);
        end
    end
    STRUC.wetstr_mc(isnan(STRUC.xhard))=nan;
    if S.struct|S.transmission|S.perm|S.revet
       %% write logfile
       % struct2log(STRUC,'STRUC','a');
    end
end

function [STRUC]=prepare_structures(S,COAST)
% function [x_hard,y_hard]=prepare_structures(S)
%
% INPUT
%     S
%         .struc
%         .x_hard
%         .y_hard
%         .XYoffset
%         .LDBstructures
%         .diffrac
%         .perm
%         .x_perm
%         .y_perm
%         .transmission  
%         .transmbwdepth 
%         .transmcrestheight 
%         .transmslope
%         .transmcrestwidth
%         .transmdir
%         .transmform
%
% OUTPUT
%     STRUC
%         .x_hard
%         .y_hard
%         .xtip
%         .ytip.
%         .hstip
%         .delta0
%         .tip_wet
%         .transmission  
%         .transmbwdepth 
%         .transmcrestheight 
%         .transmslope
%         .transmcrestwidth
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

    fprintf('  Prepare structures \n');
    STRUC=struct;
    
    % Wave blocking structure
    STRUC.x_hard=[];
    STRUC.y_hard=[];
    STRUC.s_hard=[];
    STRUC.n_hard=0;
    STRUC.type=S.structtype;
    STRUC.transmission = S.transmission;
    STRUC.transmform = S.transmform;
    STRUC.transmdir=S.transmdir;
    
    % Wave transmitting structures
    STRUC.perm=S.perm;                                                                  % switch for using hard structures
    STRUC.LDBpermeable=S.LDBpermeable;                                                              % LDB with perm structures ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
    STRUC.x_perm=[];     % for permeable structures
    STRUC.y_perm=[];     % for permeable structures
    STRUC.wavetransm=S.wavetransm;
    
    % Revetments     
    STRUC.revet=S.revet;
    STRUC.x_revet=[];     % for permeable structures
    STRUC.y_revet=[];     % for permeable structures
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
    STRUC.bypassdistribution_power=S.bypassdistribution_power;
    
    %% Coastal structures (blocking waves)
    if S.struct 
        %  add groynes/breakwaters interactively with the graphical user interface
        if strcmpi(S.LDBstructures,'manual') || strcmpi(S.LDBstructures,'interactive') 
            figure(11);
            axis equal;
            xl=xlim;yl=ylim;
            htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add structure (LMB); Next structure (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
            [x_hard,y_hard]=select_multi_polygon('k');
            set(htxt2,'Visible','off');
            STRUC.x_hard=x_hard;
            STRUC.y_hard=y_hard;

        %  try to load a file with groyne/breakwater 
        elseif ~isempty(S.LDBstructures)
            xy_hard=load(S.LDBstructures);
            STRUC.x_hard=xy_hard(:,1)'-S.XYoffset(1);
            STRUC.y_hard=xy_hard(:,2)'-S.XYoffset(2);

        % add a groyne by directly specifying the x_hard and y_hard
        elseif ~isempty(S.x_hard) && ~isempty(S.y_hard)
            STRUC.x_hard=S.x_hard(:)'-S.XYoffset(1);
            STRUC.y_hard=S.y_hard(:)'-S.XYoffset(2);
        % no groynes or breakwaters
        else
            S.struct=0;
            STRUC.x_hard=[];
            STRUC.y_hard=[];
        end
    end
    if ~isempty(STRUC.x_hard)
        STRUC.n_hard=length(find(isnan(STRUC.x_hard)))+1;
    end
    
    if STRUC.transmission == 1
        if ~isempty(S.transmfile)
            charac=load(S.transmfile);
            if size(charac(1,:),2)<STRUC.n_hard
                STRUC.transmbwdepth(1:STRUC.n_hard)=charac(1,:);
                STRUC.transmcrestheight(1:STRUC.n_hard)=charac(2,:);
                STRUC.transmslope(1:STRUC.n_hard)=charac(3,:);
                STRUC.transmcrestwidth(1:STRUC.n_hard)=charac(4,:);
                if length(charac(:,1)) == 4
                    STRUC.transmd50(1:STRUC.n_hard)=S.transmd50; % default value
                else
                    STRUC.transmd50(1:STRUC.n_hard)=charac(5,:);
                end
            elseif size(charac(1,:),2)==STRUC.n_hard
                STRUC.transmbwdepth=charac(1,:);                           
                STRUC.transmcrestheight=charac(2,:);
                STRUC.transmslope=charac(3,:);
                STRUC.transmcrestwidth=charac(4,:);
                STRUC.transmd50=charac(5,:);
            end
        elseif ~isempty(S.x_hard) && ~isempty(S.y_hard)
            if size(S.transmbwdepth,2)<STRUC.n_hard
                STRUC.transmbwdepth(1:STRUC.n_hard)=S.transmbwdepth;
            elseif size(S.transmbwdepth,2)==STRUC.n_hard
               STRUC.transmbwdepth=S.transmbwdepth;
            end
            if size(S.transmcrestheight,2)<STRUC.n_hard
                STRUC.transmcrestheight(1:STRUC.n_hard)=S.transmcrestheight;
            elseif size(S.transmcrestheight,2)==STRUC.n_hard
               STRUC.transmcrestheight=S.transmcrestheight;
            end
            if size(S.transmslope,2)<STRUC.n_hard
               STRUC.transmslope(1:STRUC.n_hard)=S.transmslope;
            elseif size(S.transmslope,2)==STRUC.n_hard
               STRUC.transmslope=S.transmslope;
            end
            if size(S.transmcrestwidth,2)<STRUC.n_hard
               STRUC.transmcrestwidth(1:STRUC.n_hard)=S.transmcrestwidth;
            elseif size(S.transmcrestwidth,2)==STRUC.n_hard
               STRUC.transmcrestwidth=S.transmcrestwidth;
            end
            if size(S.transmd50,2)<STRUC.n_hard
                STRUC.transmd50(1:STRUC.n_hard)=S.transmd50;
            elseif size(S.transmd50,2)==STRUC.n_hard
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
            STRUC.transmform=repmat(STRUC.transmform,[1,STRUC.n_hard]);
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
        if strcmpi(S.LDBpermeable,'manual') || strcmpi(S.LDBpermeable,'interactive') 
            figure(11);
            axis equal;
            xl=xlim;yl=ylim;
            htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add permeable structure (LMB); Next structure (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
            [x_perm,y_perm]=select_multi_polygon('k');
            set(htxt2,'Visible','off');
            STRUC.x_perm=x_perm;
            STRUC.y_perm=y_perm;
        
        %  try to load a file with permeable structures 
        elseif ~isempty(S.LDBpermeable)
            xy_perm=load(S.LDBpermeable);
            STRUC.x_perm=xy_perm(:,1)'-S.XYoffset(1);
            STRUC.y_perm=xy_perm(:,2)'-S.XYoffset(2);

        % add a permeable structure by directly specifying the x_perm and y_perm
        elseif ~isempty(S.x_perm) && ~isempty(S.y_perm)
            STRUC.x_perm=S.x_perm(:)'-S.XYoffset(1);
            STRUC.y_perm=S.y_perm(:)'-S.XYoffset(2);

        % no permeable structures
        else
            S.perm=0;
            STRUC.x_perm=[];
            STRUC.y_perm=[];
        end
    end
    
    %% Revetments
    STRUC.n_revet=0;
    if S.revet 
        %  add revetments interactively with the graphical user interface
        if strcmpi(S.LDBrevetments,'manual') || strcmpi(S.LDBrevetments,'interactive') 
            figure(11);
            axis equal;
            xl=xlim;yl=ylim;
            htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add revetment (LMB); Next revetment (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
            [x_revet,y_revet]=select_multi_polygon('k');
            set(htxt2,'Visible','off');
            STRUC.x_revet=x_revet;
            STRUC.y_revet=y_revet;

        %  try to load a file with revetment 
        elseif ~isempty(S.LDBrevetments)
            xy_revet=load(S.LDBrevetments);
            STRUC.x_revet=xy_revet(:,1)'-S.XYoffset(1);
            STRUC.y_revet=xy_revet(:,2)'-S.XYoffset(2);

        % add a revetment by directly specifying the x_revet and y_revet
        elseif ~isempty(S.x_revet) && ~isempty(S.y_revet)
            STRUC.x_revet=S.x_revet(:)'-S.XYoffset(1);
            STRUC.y_revet=S.y_revet(:)'-S.XYoffset(2);

        % no revetments
        else
            S.revet=0;
            STRUC.x_revet=[];
            STRUC.y_revet=[];
        end
    end
    if ~isempty(STRUC.x_revet) 
        STRUC.n_revet=sum(isnan(STRUC.x_revet))+1;
        STRUC.x_hard=[STRUC.x_hard,nan,STRUC.x_revet];
        STRUC.y_hard=[STRUC.y_hard,nan,STRUC.y_revet];
    end
    
    % Check whether to include diffraction
    STRUC.diffraction=S.diffraction;
    if isempty(STRUC.x_hard) && isempty(STRUC.x_revet) && isempty(S.x_sedlim) && isempty(S.LDBsedlim) 
        STRUC.diffraction=0; 
    end
    
    % tidy up the structures without nans at the start and end
    idnotnan=find(~isnan(STRUC.x_hard));
    if ~isempty(STRUC.x_hard)
        STRUC.x_hard=STRUC.x_hard(idnotnan(1):idnotnan(end));
        STRUC.y_hard=STRUC.y_hard(idnotnan(1):idnotnan(end));
        idnan=find(isnan(STRUC.x_hard));
        iduse=setdiff([1:length(STRUC.x_hard)],idnan(diff(idnan)==1));
        STRUC.x_hard=STRUC.x_hard(iduse);
        STRUC.y_hard=STRUC.y_hard(iduse);
    end
    
    % tidy up the permeable structures without nans at the start and end
    idnotnan=find(~isnan(STRUC.x_perm));
    if ~isempty(STRUC.x_perm)
        STRUC.x_perm=STRUC.x_perm(idnotnan(1):idnotnan(end));
        STRUC.y_perm=STRUC.y_perm(idnotnan(1):idnotnan(end));
        idnan=find(isnan(STRUC.x_perm));
        iduse=setdiff([1:length(STRUC.x_perm)],idnan(diff(idnan)==1));
        STRUC.x_perm=STRUC.x_perm(iduse);
        STRUC.y_perm=STRUC.y_perm(iduse);
    end

    % tidy up the revetments without nans at the start and end
    idnotnan=find(~isnan(STRUC.x_revet));
    if ~isempty(STRUC.x_revet)
        STRUC.x_revet=STRUC.x_revet(idnotnan(1):idnotnan(end));
        STRUC.y_revet=STRUC.y_revet(idnotnan(1):idnotnan(end));
        idnan=find(isnan(STRUC.x_revet));
        iduse=setdiff([1:length(STRUC.x_revet)],idnan(diff(idnan)==1));
        STRUC.x_revet=STRUC.x_revet(iduse);
        STRUC.y_revet=STRUC.y_revet(iduse);
    end
       
    %% IDENTIFY THE STRUCTURES WHICH ARE LOCATED IN THE SEA (OR AT LAND)
    nmc=length(COAST.x_mc)-1;
    for ist=1:length(STRUC.x_hard)
        [~,icl]=min(hypot(COAST.x_mc-STRUC.x_hard(ist),COAST.y_mc-STRUC.y_hard(ist)));
        if ~isnan(STRUC.x_hard(ist))
            im1=max(icl-1,1);
            ip1=min(icl+1,nmc+1);
            if isnan(COAST.x_mc(im1)); im1=im1-1;if im1==0;im1=2;end; end
            if isnan(COAST.x_mc(ip1)); ip1=ip1+1;if ip1>nmc+1;ip1=nmc;end; end
            dirm=360-atan2d(COAST.y_mc(ip1)-COAST.y_mc(im1),COAST.x_mc(ip1)-COAST.x_mc(im1)); 
            dirstr=atan2d(STRUC.x_hard(ist)-COAST.x_mc(icl),STRUC.y_hard(ist)-COAST.y_mc(icl));
            STRUC.wetstr_mc(ist)=int8(cosd(dirstr-dirm)>0);  % for octave, logicals cannot be assigned NaNs later on
        else
            STRUC.wetstr_mc(ist)=int8(0);
        end
    end
    STRUC.wetstr_mc(isnan(STRUC.x_hard))=nan;
    if S.struct|S.transmission|S.perm|S.revet
       %% write logfile
       % struct2log(STRUC,'STRUC','a');
    end
end

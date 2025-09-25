function [O,P]=save_shorelines(O,P,TIME,COAST,WAVE,TRANSP,STRUC,NOUR,FNOUR,GROYNE,DUNE,MUD)
% function [O,P]=save_shorelines(O,P,TIME,COAST,WAVE,TRANSP,STRUC,NOUR,FNOUR,GROYNE,DUNE,MUD)
%
% INPUT: 
%     O                    : Data structure with raw instanteneous data (where data is added to)
%     P                    : Data structure with projected and time-averaged data on a grid (where data is added to)
%         .xg              : x-coordinates of output grid on which data is projected
%         .yg              : y-coordinates of output grid on which data is projected
%         .itout           : current output timestep 
%         .cntr            : number of timesteps that have evolved for the current output timestep 'ind'
%     TIME
%         .it              : timestep index (starts at it=0 at model start)
%         .dt              : timestep [years]
%         .tc              : ratio of adaptive time step (tc=0 means using a fixed timestep)
%         .nt              : number of timesteps
%         .tnow            : current model time [days in datenum format]
%         .timenum0        : model start time [days in datenum format]
%         .itout           : timestep index of next output moment   
%         .storageinterval : interval between output moments [days]
%     COAST
%         .x_mc            : x-coordinates of the coastline of all elements
%         .y_mc            : y-coordinates of the coastline of all elements
%         .x_mc1           : x-coordinates of the initial coastline of all elements at t0
%         .y_mc1           : y-coordinates of the initial coastline of all elements at t0
%         .dSds            : coastline change rate 
%         .PHIc            : coastline orientation of the current coastline [°N]
%         .QS              : transport along the coast [m3/yr]
%         .QSmax           : maximum theoretical transport capacity along the coast [m3/yr]
%     WAVE
%         .HSo             : offshore wave height along the coast [m]
%         .PHIo            : offshore wave direction along the coast [°N]
%         .tper            : wave period along the coast [s]
%         .PHIf            : lower-shoreface orientation [°N]
%         .HS              : significant wave height at depth-of-closure [m]
%         .PHI             : wave direction at depth-of-closure [°N] 
%         .dPHI            : relative angle of nearshore waves at depth-of-closure w.r.t. coastline [°]
%         .HSbr            : significant wave height at point of breaking [m]
%         .PHIbr           : wave direction at point of breaking [°N] 
%         .dPHIbr          : relative angle of waves at point of breaking w.r.t. coastline [°]  
%         .hbr             : waves at point of breaking
%     STRUC                : structures
%         .xhard           : x-coordinates of the hard structures
%         .yhard           : y-coordinates of the hard structures
%     NOUR
%         .xnour           : x-coordinates of the nourishment locations
%         .ynour           : y-coordinates of the nourishment locations
%     GROYNE               : groyne locations
%         .x               : x-coordinates of the groynes
%         .y               : y-coordinates of the groynes
% OUTPUT:
%     O                    : Data-structure with raw instanteneous data, with corresponding fields to the input variables
%         .it              : timestep index (starts at it=0 at model start)
%         .dt              : timestep [years]
%         .tc              : ratio of adaptive time step (tc=0 means using a fixed timestep)
%         .nt              : number of timesteps
%         .time            : current model time [days w.r.t. start of model]
%         .timenum         : current model time [days in datenum format]
%         .n               : number of cells of all coastline elements [-]
%         .x               : x-coordinates of all coastline elements [m], where elements are separated by NaNs
%         .y               : y-coordinates of all coastline elements [m], where elements are separated by NaNs
%         .h0              : active height of the coastal grid cells
%         .distx           : alongshore distance along the coastline elements [m]
%         .PHIc            : coastline orientation of the current coastline at transport points [°N]
%         .PHIcxy          : coastline orientation of the current coastline [°N]
%         .PHIf            : lower-shoreface orientation [°N]
%         .distQS          : alongshore transport along the coastline elements [m]
%         .QS              : transport along the coast [m3/yr]
%         .TP              : wave period along the coast [s]
%         .HSo             : offshore wave height along the coast [m]
%         .HStdp           : significant wave height at depth-of-closure [m]
%         .HSbr            : significant wave height at point of breaking [m]
%         .PHIo            : wave direction along the coast at offshore point [°N]
%         .PHItdp          : wave direction at depth-of-closure [°N] 
%         .PHIbr           : wave direction at point of breaking [°N] 
%         .wberm           : width of the berm/beach [m]
%         .qs              : dune erosion volume change that increases the beach width [m3/m/yr]
%         .ql              : dune erosion volume change that does not increase the beach width [m3/m/yr]
%         .qw              : wind transport from the beach to the dune [m3/m/yr]
%         .R               : runup level for dune erosion [m]
%         .SWL             : still water level for dune erosion [m]
%         .PHIcxy          : shore-normal orientation at coastline points [°N]
%         .xdune           : x-coordinates of the dunes [m]
%         .ydune           : y-coordinates of the dunes [m]
%         .dfelev          : duen foot height [m MSL]
%         .dcelev          : dune crest height [m MSL]
%         .hd0             : active height of the dune grid cells [m]
%         .Bf              : mud flat width [m]
%         .Bm              : mangrove width [m]
%         .Bfm             : colonizing mangrove width [m]
%         .xhard           : x-coordinates of the hard structures
%         .yhard           : y-coordinates of the hard structures
%         .nhard           : number of structures
%         .xnour           : x-coordinates of the nourishment locations
%         .ynour           : y-coordinates of the nourishment locations
%         .nnour           : number of nourishments
%         .x_fnour         : x-coordinates of shoreface nourishment locations
%         .y_fnour         : y-coordinates of shoreface nourishment locations
%         .n_fnour         : number of shoreface nourishments
%         .V_fnour_t       : volume of shoreface nourishments
%         .q_fnour_t       : discharge rate of shoreface nourishments
%     P                    : Data structure with projected on a grid which is for the waves and transports also time-averaged data in-between output timesteps, with the following fields:
%         .it              : timestep index (starts at it=0 at model start)
%         .itout           : index for the next output step of the P-structure
%         .cntr            : reset counter for moving average of current output step of the P-structure
%         .tc              : ratio of adaptive time step (tc=0 means using a fixed timestep)
%         .dt              : timestep [years]
%         .nt              : number of timesteps
%         .time            : model time at output timesteps [days w.r.t. start of model]
%         .timenum         : model time at output timesteps [days in datenum format]
%         .xg              : x-coordinates of the output grids [m], where elements are separated by NaNs
%         .yg              : y-coordinates of the output grids [m], where elements are separated by NaNs
%         .xc              : x-coordinates of the fitted coastline at each output timestep w.r.t. the output grid [m]
%         .yc              : y-coordinates of the fitted coastline at each output timestep w.r.t. the output grid [m]
%         .zc              : cross-shore position of the fitted coastline at each output timestep w.r.t. shorenormal transects of the output grid [m]
%         .volc            : cross-shore volume of the fitted coastline at each output timestep w.r.t. shorenormal transects of the output grid [m]
%         .h0              : active height of the coastal grid cells [m]
%         .xd              : x-coordinates of the fitted dune-line at each output timestep w.r.t. the output grid [m]
%         .yd              : y-coordinates of the fitted dune-line at each output timestep w.r.t. the output grid [m]
%         .zd              : cross-shore position of the fitted dune-line at each output timestep w.r.t. shorenormal transects of the output grid [m]
%         .vold            : cross-shore volume of the fitted dune-line at each output timestep w.r.t. shorenormal transects of the output grid [m]
%         .hd0             : active height of the dune grid cells [m]
%         .PHIc            : orientation of the coastline at transport points along the output grid [°N]
%         .PHIcxy          : orientation of the coastline along the output grid [°N]
%         .PHIf            : lower-shoreface orientation along the output grid [°N]
%         .TP              : wave period along the output grid [s]
%         .HSo             : offshore wave height along the output grid [m]
%         .HStdp           : significant wave height at depth-of-closure along the output grid [m]
%         .HSbr            : significant wave height at point of breaking along the output grid [m]
%         .PHIo            : wave direction at offshore points along the output grid [°N] 
%         .PHItdp          : wave direction at depth-of-closure along the output grid [°N] 
%         .PHIbr           : wave direction at point of breaking along the output grid [°N] 
%         .QS              : transport along the coast [m3/yr]
%         .QSmax           : maximum theoretical transport capacity along the coast [m3/yr]
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

    % only save data when a stoarge time-interval has passed
    if isempty(O.storagedate)
        storagedatenum=unique([TIME.timenum0:O.storageinterval:TIME.tend,TIME.tend]);
    else
        storagedatenum=unique(datenum(O.storagedate));
    end
    storedata=0;
    if TIME.dtsteps==0 && O.itout<=length(storagedatenum)
        storedata=TIME.tnow>=storagedatenum(O.itout);
    end
    if (TIME.tnow+TIME.dt)>=TIME.tend 
        storedata=1;
    end
    
    %% export directly the data to O structure (with variable grid size)
    if storedata
        it=TIME.it;
        dt=TIME.dt;
        tc=TIME.tc;
        nt=TIME.nt;
        x_mc=COAST.x1_mc;
        y_mc=COAST.y1_mc;
        PHIc_mc=COAST.PHIc1_mc;
        PHIcxy_mc=COAST.PHIcxy1_mc;
        if (TIME.tnow+TIME.dt)>=TIME.tend
            x_mc=COAST.x_mc;
            y_mc=COAST.y_mc;
            PHIc_mc=COAST.PHIc_mc;
            PHIcxy_mc=COAST.PHIcxy_mc;
        end
        h0_mc=COAST.h0_mc;
        QS=TRANSP.QS_mc;                    
        HSo=WAVE.HSo_mc;
        PHIo=WAVE.PHIo_mc;
        TP=WAVE.TP_mc;
        PHIf=COAST.PHIf_mc;                 % offshore waves & lower-shoreface orientation
        HStdp=WAVE.HStdp_mc;
        PHItdp=WAVE.PHItdp_mc;              % nearshore waves at depth-of-closure
        HSbr=WAVE.HSbr_mc;
        PHIbr=WAVE.PHIbr_mc;
        %hbr=WAVE.hbr_mc;                   % depth at point of breaking
        if MUD.used
            Bf_mc=COAST.Bf_mc;
            Bm_mc=COAST.Bm_mc;
            Bfm_mc=COAST.Bfm_mc;
        end
        if DUNE.used
            wberm_mc=COAST.wberm_mc;
            qs_mc=COAST.qs_mc;
            ql_mc=COAST.ql_mc;
            qw_mc=COAST.qw_mc;
            R_mc=COAST.R_mc;
            SWL_mc=COAST.SWL_mc;
            xdune_mc = COAST.xdune_mc;
            ydune_mc = COAST.ydune_mc;
            dfelev_mc = COAST.dfelev_mc;                  % dune foot height [m MSL]
            dcelev_mc = COAST.dcelev_mc;                  % dune crest height [m MSL]
            hd0_mc = dcelev_mc-dfelev_mc;                 % active height of the dune grid cells [m]
        end
        xhard=STRUC.xhard;
        yhard=STRUC.yhard;                  % structures
        xnour=NOUR.xnour;
        ynour=NOUR.ynour;                   % nourishment locations
        x_groyne=GROYNE.x;
        y_groyne=GROYNE.y;                  % groyne locations on shoreline
        if FNOUR.fnourish == 1
            x_fnour=FNOUR.x;                % shoreface nourishment info 
            y_fnour=FNOUR.y;
            V_fnour=FNOUR.Vt(:,end); 
            q_fnour=FNOUR.q_tot_mc; 
        end
        
        %% compute distance x,y and transport locations
        dx=x_mc(2:end)-x_mc(1:end-1);dx(isnan(dx))=0;
        dy=y_mc(2:end)-y_mc(1:end-1);dy(isnan(dy))=0;
        distx = [0,cumsum((dx.^2+dy.^2).^0.5)];
        distQS = (distx(1:end-1)+distx(2:end))/2;
        
        % store data in structure (which is automatically resized if needed)
        [O.it]        = addVECTOR(O.it, it);
        [O.dt]        = addVECTOR(O.dt, dt);
        [O.tc]        = addVECTOR(O.dt, tc);
        [O.nt]        = addVECTOR(O.nt, nt);
        [O.time]      = addVECTOR(O.time, TIME.tnow-TIME.timenum0);
        [O.timenum]   = addVECTOR(O.timenum, TIME.tnow);
        
        % parameters after coastline change
        [O.x]         = addVECTOR(O.x, x_mc');
        [O.y]         = addVECTOR(O.y, y_mc');
        [O.distx]     = addVECTOR(O.distx, distx');
        [O.h0]        = addVECTOR(O.h0, h0_mc'); 

        % coastline orientation and active height
        [O.PHIc]      = addVECTOR(O.PHIc, PHIc_mc');
        [O.PHIcxy]    = addVECTOR(O.PHIcxy, PHIcxy_mc'); 
        [O.PHIf]      = addVECTOR(O.PHIf, PHIf);

        % transport parameters
        [O.distQS]    = addVECTOR(O.distQS, distQS');
        [O.QS]        = addVECTOR(O.QS, QS');

        % offshore waves
        [O.TP]        = addVECTOR(O.TP, TP);
        [O.HSo]       = addVECTOR(O.HSo, HSo);
        [O.PHIo]      = addVECTOR(O.PHIo, PHIo);

        % nearshore waves at depth-of-closure
        [O.HStdp]     = addVECTOR(O.HStdp, HStdp');
        [O.PHItdp]    = addVECTOR(O.PHItdp, PHItdp');
        
        % waves at point of breaking
        [O.HSbr]      = addVECTOR(O.HSbr, HSbr);
        [O.PHIbr]     = addVECTOR(O.PHIbr, PHIbr);
        %[O.hbr]       = addVECTOR(O.hbr, hbr);

        if DUNE.used
            [O.wberm] = addVECTOR(O.wberm, wberm_mc'); 
            [O.qs]    = addVECTOR(O.qs, qs_mc'); 
            [O.ql]    = addVECTOR(O.ql, ql_mc'); 
            [O.qw]    = addVECTOR(O.qw, qw_mc'); 
            [O.R]     = addVECTOR(O.R, R_mc'); 
            [O.SWL]   = addVECTOR(O.SWL, SWL_mc'); 
            [O.xdune] = addVECTOR(O.xdune, xdune_mc'); 
            [O.ydune] = addVECTOR(O.ydune, ydune_mc'); 
            [O.dfelev] = addVECTOR(O.dfelev, dfelev_mc');                  % dune foot height [m MSL]
            [O.dcelev] = addVECTOR(O.dcelev, dcelev_mc');                  % dune crest height [m MSL]
            [O.hd0] = addVECTOR(O.hd0, hd0_mc');                           % active height of the dune grid cells [m]
        end
        if MUD.used
            [O.Bf]    = addVECTOR(O.Bf, Bf_mc');
            [O.Bm]    = addVECTOR(O.Bm, Bm_mc');
            [O.Bfm]   = addVECTOR(O.Bfm, Bfm_mc');
        end

        % structures
        [O.xhard]     = addVECTOR(O.xhard, xhard');
        [O.yhard]     = addVECTOR(O.yhard, yhard');
        [O.nhard]     = addVECTOR(O.nhard, length(yhard));
        
        % nourishments
        [O.xnour]    = addVECTOR(O.xnour, xnour');
        [O.ynour]    = addVECTOR(O.ynour, ynour');
        [O.nnour]    = addVECTOR(O.nnour, length(ynour));
        
        % groyne locations at shoreline
        [O.x_groyne]  = addVECTOR(O.x_groyne, x_groyne(:));
        [O.y_groyne]  = addVECTOR(O.y_groyne, y_groyne(:));
        
        % shoreface nourishments (added by Anne)
        if FNOUR.fnourish == 1
            [O.x_fnour]    = addVECTOR(O.x_fnour, x_fnour(:));
            [O.y_fnour]    = addVECTOR(O.y_fnour, y_fnour(:));
            [O.n_fnour]    = FNOUR.n;
            [O.V_fnour_t]  = addVECTOR(O.V_fnour_t, V_fnour(:)); % remaining volume after each time step [m3]
            [O.q_fnour_t]  = addVECTOR(O.q_fnour_t, q_fnour(:)); % [m3/m/yr]
        end 
        
        O.itout=O.itout+1;
    end
    
    %% Project output data on a grid to P structure (with fixed grid position & mean of transport and waves)
    outputgridnr=0;
    if isfield(P,'itout')
        outputgridnr=length(P);
    end
    if TIME.dtsteps==0
        for pp=1:outputgridnr
            % current output timestep
            ind=P(pp).itout;
            % number of timesteps that have evolved for the current output timestep 'ind'
            cntr=P(pp).cntr;
            
            % OUTPUT: projection grid used for exporting data
            xp=P(pp).xg;
            yp=P(pp).yg;
            if ~isempty(P(pp).xc)
                xc=P(pp).xc(:,end);
                yc=P(pp).yc(:,end);
            else
                xc=xp(:);
                yc=yp(:);
            end
            
            % xq grid stored after the collect variables step
            xw=COAST.xq1_mc; 
            xw2=COAST.x1_mc; 
            yw=COAST.yq1_mc;
            yw2=COAST.y1_mc;
            % interpolate data to the projection grid (xp,yp)
            fields1={'TP','HSo','HStdp','HSbr','QS','QSmax','h0'};
            fields2={'PHIc','PHIf','PHIo','PHItdp','PHIbr','PHIcxy'};
            if DUNE.used
                fields1=[fields1,{'dfelev','dcelev','wberm','qs','ql','qw','R','SWL'}];
            end
            if MUD.used
                fields1=[fields1,{'Bf','Bm','Bfm'}];
            end
            nmax=length(fields1);
            
            [scalars,vectors]=get_gridprojection(COAST,WAVE,TRANSP,xc,yc,xw,yw,fields1(1:6),fields2(1:5),'weighted_distance');
            [scalars2,vectors2]=get_gridprojection(COAST,WAVE,TRANSP,xc,yc,xw2,yw2,fields1(7:nmax),fields2(6),'weighted_distance');
            for ff=1:length(fields1)
                if ~isfield(scalars,fields1{ff})
                    scalars.(fields1{ff})=scalars2.(fields1{ff});
                end
            end
            for ff=1:length(fields2)
                if ~isfield(vectors,fields2{ff})
                    vectors.(fields2{ff})=vectors2.(fields2{ff});
                end
            end
            
            % initialize a new column of data if cntr=1
            if cntr==1
                for ff=1:length(fields1)
                P(pp).(fields1{ff})(:,ind)=zeros(length(P(pp).xg),1);
                end
                for ff=1:length(fields2)
                P(pp).(fields2{ff})(:,ind)=zeros(length(P(pp).xg),1);
                end
            end
            
            % continuously update the data stored as a moving average in the P structure        
            wghtNR=[cntr-1,1]/cntr;
            HStdp0=max(P(pp).HStdp(:,ind),1e-6);
            for ff=1:length(fields1)
                P(pp).(fields1{ff})(:,ind) = P(pp).(fields1{ff})(:,ind)*wghtNR(1) + scalars.(fields1{ff})(:)*wghtNR(2);
            end
            wghtHS=[HStdp0.^2,scalars.HStdp(:).^2];
            for ff=1:length(fields2)
                sdr = sind(P(pp).(fields2{ff})(:,ind)).*wghtNR(1).*wghtHS(:,1) + sind(vectors.(fields2{ff})(:))*wghtNR(2).*wghtHS(:,2);
                cdr = cosd(P(pp).(fields2{ff})(:,ind)).*wghtNR(1).*wghtHS(:,1) + cosd(vectors.(fields2{ff})(:))*wghtNR(2).*wghtHS(:,2);
                P(pp).(fields2{ff})(:,ind) = mod(atan2d(sdr,cdr),360);
            end
            
            % store time variables
            P(pp).tc=TIME.tc;
            P(pp).it(1,ind)=TIME.it;
            P(pp).dt(1,ind)=TIME.dt;
            P(pp).nt(1,ind)=TIME.nt;
            P(pp).time(1,ind)=TIME.tnow-TIME.timenum0;
            P(pp).timenum(1,ind)=TIME.tnow;
            if ind>1 && TIME.tc~=0
            P(pp).dt(1,ind)=(TIME.tnow-P(pp).timenum(1,ind-1))/P(pp).cntr/365;
            end
            
            % store the cross-shore location of the coastline w.r.t. 
            if storedata
                Lcrit=5000;
                x_mc=interpNANs(COAST.x1_mc);
                y_mc=interpNANs(COAST.y1_mc);
                
                % extend coastline with 1% of a grid cell length
                if ~COAST.cyclic && COAST.n_mc==1
                    dx=diff(x_mc);
                    dy=diff(y_mc);
                    x_mc=[x_mc(1)-dx(1)/100,x_mc,x_mc(end)+dx(end)/100];
                    y_mc=[y_mc(1)-dy(1)/100,y_mc,y_mc(end)+dy(end)/100];
                end
                
                % add coastline position w.r.t. reference line (zc) and as xy-coordinates (xc,yc)
                [cmin,xcr,ycr]=get_polydistance(xp,yp,x_mc,y_mc,Lcrit);
                P(pp).zc=[P(pp).zc,cmin(:)];
                P(pp).xc=[P(pp).xc,xcr(:)];
                P(pp).yc=[P(pp).yc,ycr(:)];
                P(pp).volc=[P(pp).volc,cmin(:).*P(pp).h0(:,end)];
                
                % add dune position w.r.t. the reference line (zd) and as xy-coordinates (xd, yd)
                if DUNE.used
                    xdune_mc=interpNANs(COAST.xdune_mc);
                    ydune_mc=interpNANs(COAST.ydune_mc);
                    
                    % extend duneline with 10% of a grid cell length
                    if ~COAST.cyclic && COAST.n_mc==1
                        dx=diff(xdune_mc);
                        dy=diff(ydune_mc);
                        xdune_mc=[xdune_mc(1)-dx(1)/10,xdune_mc,xdune_mc(end)+dx(end)/10];
                        ydune_mc=[ydune_mc(1)-dy(1)/10,ydune_mc,ydune_mc(end)+dy(end)/10];
                    end

                    % add duneline position w.r.t. reference line (zd) and as xy-coordinates (xd,yd)
                    [dmin,xdr,ydr]=get_polydistance(xp,yp,xdune_mc,ydune_mc,Lcrit);
                    P(pp).zd=[P(pp).zd,dmin(:)];
                    P(pp).xd=[P(pp).xd,xdr(:)];
                    P(pp).yd=[P(pp).yd,ydr(:)];
                    P(pp).hd0=[P(pp).hd0,P(pp).dcelev(:,end)-P(pp).dfelev(:,end)];
                    P(pp).vold=[P(pp).vold,dmin(:).*P(pp).hd0(:,end)];
                end
                
                % reset counter for moving average of current output step of the P-structure
                P(pp).cntr=1;
                % set index for the next output step of the P-structure
                P(pp).itout=P(pp).itout+1;
            else
                % count the number of times or moving average of current output step of the P-structure
                P(pp).cntr=P(pp).cntr+1;
            end
        end
    end
    
    %% Store data at predefined profiles
    % these data can be plotted over time
    if ~isempty(O.xyprofiles)
        xp=O.xyprofiles(:,1);
        yp=O.xyprofiles(:,2);
        
        % idp=[];
        % for pp=1:length(xp)
        %     dist=((x_mc-xp(pp)).^2+(y_mc-yp(pp)).^2).^0.5;
        %     idp(pp)=find(dist==min(dist),1);
        % end
        
        % interpolate coastline data to the projection grid (xp,yp)      
        Lcrit=2000;
        [cmin,xcr,ycr]=get_polydistance(xp,yp,COAST.x_mc,COAST.y_mc,Lcrit);
        [O.timenum_profile] = addVECTOR(O.timenum_profile,TIME.tnow);
        [O.adt_profile] = addVECTOR(O.adt_profile, TIME.dt);
        [O.it_profile] = addVECTOR(O.it_profile, TIME.it);
        [O.c_profile] = addVECTOR(O.c_profile, cmin(:));
        [O.xc_profile] = addVECTOR(O.xc_profile, xcr(:));
        [O.yc_profile] = addVECTOR(O.yc_profile, ycr(:));
        if DUNE.used
            % interpolate dune data to the projection grid (xp,yp)      
            [dmin,xdr,ydr]=get_polydistance(xp,yp,COAST.xdune_mc,COAST.ydune_mc,Lcrit);
            [O.d_profile] = addVECTOR(O.d_profile, dmin(:));
            [O.xd_profile] = addVECTOR(O.xd_profile, xdr(:));
            [O.yd_profile] = addVECTOR(O.yd_profile, ydr(:));
        end
    end
    
    %% Store all shoreline data to file
    if storedata
        if ~isoctave, warning off, end
        save(fullfile(pwd,O.outputdir,O.outputfile),'O','P');
        if ~isoctave, warning on, end
    end
end

% function to extend the size of a matrix
function [c]=addVECTOR(a,b)
b=b(:);
S1 = size(a);
S2 = length(b);
if S2>S1(1)
   a = [a;nan(S2-S1(1),S1(2))];
elseif S1(1)>S2
   b = [b;nan(S1(1)-S2,1)];
end
c = [a,b];
end
function [O,P,itout]=save_shorelines(O,P,S,TIME,COAST,WAVE,TRANSP,STRUC,NOUR,FNOUR,GROYNE,DUNE,MUD)
% function [O,P,itout]=save_shorelines(O,P,S,TIME,COAST,WAVE,TRANSP,STRUC,NOUR,FNOUR,GROYNE,DUNE,MUD)
%
% INPUT: 
%     O                    : Data structure with raw instanteneous data (where data is added to)
%     P                    : Data structure with projected and time-averaged data on a grid (where data is added to)
%         .xg              : x-coordinates of output grid on which data is projected
%         .yg              : y-coordinates of output grid on which data is projected
%         .itout           : current output timestep 
%         .cntr            : number of timesteps that have evolved for the current output timestep 'ind'
%     S                    : Model input data structure
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
%         .timenum         : current model time [days in datenum format]
%         .n               : number of cells of all coastline elements [-]
%         .x               : x-coordinates of all coastline elements [m], where elements are separated by NaNs
%         .y               : y-coordinates of all coastline elements [m], where elements are separated by NaNs
%         .wberm           : width of the berm/beach [m]
%         .qs              : dune erosion volume change that increases the beach width [m3/m/yr]
%         .ql              : dune erosion volume change that does not increase the beach width [m3/m/yr]
%         .qw              : wind transport from the beach to the dune [m3/m/yr]
%         .R               : runup level for dune erosion [m]
%         .SWL             : still water level for dune erosion [m]
%         .PHIcxy          : shore-normal orientation at coastline points [°N]
%         .xdune           : x-coordinates of the dunes [m]
%         .ydune           : y-coordinates of the dunes [m]
%         .Bf              : mud flat width [m]
%         .Bm              : mangrove width [m]
%         .Bfm             : colonizing mangrove width [m]
%         .distx           : alongshore distance along the coastline elements [m]
%         .distQS          : alongshore transport along the coastline elements [m]
%         .dSds            : coastline change rate
%         .n1              : number of cells of all coastline elements at t0 [-]
%         .x1              : x-coordinates of all coastline elements at t0 [m]
%         .y1              : y-coordinates of all coastline elements at t0 [m]
%         .distx1          : alongshore distance along the coastline elements at t0 [m]
%         .distQS1         : alongshore transport along the coastline elements at t0 [m]
%         .PHIc            : coastline orientation of the current coastline [°N]
%         .QS              : transport along the coast [m3/yr]
%         .HSo             : offshore wave height along the coast [m]
%         .PHIo            : offshore wave direction along the coast [°N]
%         .TP              : wave period along the coast [s]
%         .PHIf            : lower-shoreface orientation [°N]
%         .HS              : significant wave height at depth-of-closure [m]
%         .PHI             : wave direction at depth-of-closure [°N] 
%         .dPHI            : relative angle of nearshore waves at depth-of-closure w.r.t. coastline [°]
%         .HSbr            : significant wave height at point of breaking [m]
%         .PHIbr           : wave direction at point of breaking [°N] 
%         .dPHIbr          : relative angle of waves at point of breaking w.r.t. coastline [°] 
%         .hbr             : depth at point of breaking [m]
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
%         .timenum         : model time at output timesteps [days in datenum format]
%         .xg              : x-coordinates of the output grids [m], where elements are separated by NaNs
%         .yg              : y-coordinates of the output grids [m], where elements are separated by NaNs
%         .zg              : cross-shore position of the fitted coastline at each output timestep w.r.t. shorenormal transects of the output grid [m]
%         .xc              : x-coordinates of the fitted coastline at each output timestep w.r.t. the output grid [m]
%         .yc              : y-coordinates of the fitted coastline at each output timestep w.r.t. the output grid [m]
%         .TP              : wave period along the output grid [s]
%         .HSo             : offshore wave height along the output grid [m]
%         .HStdp           : significant wave height at depth-of-closure along the output grid [m]
%         .HSbr            : significant wave height at point of breaking along the output grid [m]
%         .PHIc            : orientation of the coastline along the output grid [°N]
%         .PHIf            : lower-shoreface orientation along the output grid [°N]
%         .PHIo            : wave direction at offshore along the output grid [°N] 
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
    timenum1=TIME.tnow;
    timenum0=TIME.timenum0;
    itout=TIME.itout;
    storageinterval=O.storageinterval;        % time parameters
    if TIME.tc==0
        TIME.adt=TIME.dt;
    end
    storedata=(timenum1-timenum0)>=itout*storageinterval || (TIME.tnow+TIME.adt)>=TIME.tend;
    
    %% export directly the data to O structure (with variable grid size)
    if storedata
        it=TIME.it;
        dt=TIME.dt;
        tc=TIME.tc;
        nt=TIME.nt;
        x_mc=COAST.x_mc;
        y_mc=COAST.y_mc;
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
            PHIcxy_mc=COAST.PHIcxy_mc;
            xdune = COAST.xdune;
            ydune = COAST.ydune;
        end
        dSds=COAST.dSds_mc;                 % x,y (after coastline update) and coastline change
        x_mc1=COAST.x1_mc;
        y_mc1=COAST.y1_mc;
        PHIc=COAST.PHIc_mc;
        QS=TRANSP.QS_mc;                    % x,y (before coastline update) and corresponding QS
        HSo=WAVE.HSo_mc;
        PHIo=WAVE.PHIo_mc;
        TP=WAVE.TP_mc;
        PHIf=COAST.PHIf_mc;                 % offshore waves & lower-shoreface orientation
        HS=WAVE.HStdp_mc;
        PHI=WAVE.PHItdp_mc;
        dPHI=WAVE.dPHItdp_mc;               % nearshore waves at depth-of-closure
        HSbr=WAVE.HSbr_mc;
        PHIbr=WAVE.PHIbr_mc;
        dPHIbr=WAVE.dPHIbr_mc;
        hbr=WAVE.hbr_mc;                    % waves at point of breaking
        xhard=STRUC.xhard;
        yhard=STRUC.yhard;                  % structures
        xnour=NOUR.xnour;
        ynour=NOUR.ynour;                 % nourishment locations
        x_groyne=GROYNE.x;
        y_groyne=GROYNE.y;                  % groyne locations on shoreline
        if FNOUR.fnourish == 1
            x_fnour=FNOUR.x;                % shoreface nourishment info 
            y_fnour=FNOUR.y;
            V_fnour=FNOUR.Vt(:,end); 
            q_fnour=FNOUR.q_tot_mc; 
        end

        %% INITIALIZE OUTPUT FIELDS
        nans=find(isnan(x_mc));
        nans1=find(isnan(x_mc1));
        n_mc=length(nans)+1;
        n_mc1=length(nans1)+1;

        %% compute distance x,y and transport locations
        dx=x_mc(2:end)-x_mc(1:end-1);dx(isnan(dx))=0;
        dy=y_mc(2:end)-y_mc(1:end-1);dy(isnan(dy))=0;
        distx = [0,cumsum((dx.^2+dy.^2).^0.5)];
        distQS = (distx(1:end-1)+distx(2:end))/2;
        dx=x_mc1(2:end)-x_mc1(1:end-1);dx(isnan(dx))=0;
        dy=y_mc1(2:end)-y_mc1(1:end-1);dy(isnan(dy))=0;
        distx1 = [0,cumsum((dx.^2+dy.^2).^0.5)];
        distQS1 = (distx1(1:end-1)+distx1(2:end))/2;

        % store data in structure (which is automatically resized if needed)
        [O.it]        = addVECTOR(O.it, it);
        [O.dt]        = addVECTOR(O.dt, dt);
        [O.tc]        = addVECTOR(O.dt, tc);
        [O.nt]        = addVECTOR(O.nt, nt);
        [O.timenum]   = addVECTOR(O.timenum, timenum1);

        % parameters after coastline change
        [O.n]         = addVECTOR(O.n, n_mc');
        [O.x]         = addVECTOR(O.x, x_mc');
        [O.y]         = addVECTOR(O.y, y_mc');
        if DUNE.used
            [O.wberm] = addVECTOR(O.wberm, wberm_mc'); 
            [O.qs]    = addVECTOR(O.qs, qs_mc'); 
            [O.ql]    = addVECTOR(O.ql, ql_mc'); 
            [O.qw]    = addVECTOR(O.qw, qw_mc'); 
            [O.R]     = addVECTOR(O.R, R_mc'); 
            [O.SWL]   = addVECTOR(O.SWL, SWL_mc'); 
            [O.PHIcxy]= addVECTOR(O.PHIcxy, PHIcxy_mc'); 
            [O.xdune] = addVECTOR(O.xdune, xdune'); 
            [O.ydune] = addVECTOR(O.ydune, ydune'); 
        end
        if MUD.used
            [O.Bf]    = addVECTOR(O.Bf, Bf_mc');
            [O.Bm]    = addVECTOR(O.Bm, Bm_mc');
            [O.Bfm]   = addVECTOR(O.Bfm, Bfm_mc');
        end
        [O.distx]     = addVECTOR(O.distx, distx');
        [O.distQS]    = addVECTOR(O.distQS, distQS');
        [O.dSds]      = addVECTOR(O.dSds, dSds');

        % parameters before coastline change
        [O.n1]        = addVECTOR(O.n1, n_mc1');
        [O.x1]        = addVECTOR(O.x1, x_mc1');
        [O.y1]        = addVECTOR(O.y1, y_mc1');
        [O.distx1]    = addVECTOR(O.distx1, distx1');
        [O.distQS1]   = addVECTOR(O.distQS1, distQS1');
        
        % coastline orientation and transport
        [O.PHIc]      = addVECTOR(O.PHIc, PHIc');
        [O.QS]        = addVECTOR(O.QS, QS');

        % offshore waves
        [O.HSo]       = addVECTOR(O.HSo, HSo);
        [O.PHIo]      = addVECTOR(O.PHIo, PHIo);
        [O.TP]        = addVECTOR(O.TP, TP);
        [O.PHIf]      = addVECTOR(O.PHIf, PHIf);

        % nearshore waves at depth-of-closure
        [O.HS]        = addVECTOR(O.HS, HS');
        [O.PHI]       = addVECTOR(O.PHI, PHI');
        [O.dPHI]      = addVECTOR(O.dPHI, dPHI');

        % waves at point of breaking
        [O.HSbr]      = addVECTOR(O.HSbr, HSbr);
        [O.PHIbr]     = addVECTOR(O.PHIbr, PHIbr);
        [O.dPHIbr]    = addVECTOR(O.dPHIbr, dPHIbr);
        [O.hbr]       = addVECTOR(O.hbr, hbr);

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

        itout=itout+1;
    end

    %% Project output data on a grid to P structure (with fixed grid position & mean of transport and waves)
    outputgridnr=0;
    if isfield(P,'itout')
        outputgridnr=length(P);
    end
    for pp=1:outputgridnr
        % current output timestep
        ind=P(pp).itout;
        % number of timesteps that have evolved for the current output timestep 'ind'
        cntr=P(pp).cntr;
        
        % OUTPUT: projection grid used for exporting data
        xp=P(pp).xg;
        yp=P(pp).yg;
        % xq grid stored after the collect variables step
        xw=COAST.xq1_mc; 
        yw=COAST.yq1_mc;
        % interpolate data to the projection grid (xp,yp)
        fields1={'TP','HSo','HStdp','HSbr','QS','QSmax'};
        fields2={'PHIc','PHIf','PHIo','PHItdp','PHIbr'};
        [scalars,vectors,idgrid]=get_gridprojection(COAST,WAVE,TRANSP,xp,yp,xw,yw,fields1,fields2,'weighted_distance');
        
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
        HStdp0=P(pp).HStdp(:,ind);
        for ff=1:length(fields1)
            P(pp).(fields1{ff})(:,ind) = P(pp).(fields1{ff})(:,ind)*wghtNR(1) + scalars.(fields1{ff})(:)*wghtNR(2);
        end
        wghtHS=[HStdp0.^2,scalars.HStdp(:).^2];
        for ff=1:length(fields2)
            sdr = sind(P(pp).(fields2{ff})(:,ind)).*wghtNR(1).*wghtHS(:,1) + sind(vectors.(fields2{ff})(:))*wghtNR(2).*wghtHS(:,2);
            cdr = cosd(P(pp).(fields2{ff})(:,ind)).*wghtNR(1).*wghtHS(:,1) + cosd(vectors.(fields2{ff})(:))*wghtNR(2).*wghtHS(:,2);
            P(pp).(fields2{ff})(:,ind) = mod(atan2d(cdr,sdr),360);
        end
        
        % store time variables
        P(pp).tc=TIME.tc;
        P(pp).it(1,ind)=TIME.it;
        P(pp).dt(1,ind)=TIME.dt;
        P(pp).nt(1,ind)=TIME.nt;
        P(pp).timenum(1,ind)=TIME.tnow;
        if ind>1 && TIME.tc~=0
        P(pp).dt(1,ind)=(TIME.tnow-P(pp).timenum(1,ind-1))/P(pp).cntr/365;
        end
        
        % store the cross-shore location of the coastline w.r.t. 
        if storedata
            Lcrit=5000;
            x_mc=COAST.x1_mc;
            y_mc=COAST.y1_mc;
            %try
            %    %extend with 1% of a grid cell
            %    dx=diff(x_mc);
            %    dy=diff(y_mc);
            %    x_mc=[x_mc(1)-dx(1)/100,x_mc,x_mc(end)+dx(end)/100];
            %    y_mc=[y_mc(1)-dy(1)/100,y_mc,y_mc(end)+dy(end)/100];
            %end
            [dmin,xcr,ycr]=get_polydistance(xp,yp,x_mc,y_mc,Lcrit);
            P(pp).zg=[P(pp).zg,dmin(:)];
            P(pp).xc=[P(pp).xc,xcr(:)];
            P(pp).yc=[P(pp).yc,ycr(:)];
            % reset counter for moving average of current output step of the P-structure
            P(pp).cntr=1;
            % set index for the next output step of the P-structure
            P(pp).itout=P(pp).itout+1;
        else
            % count the number of times or moving average of current output step of the P-structure
            P(pp).cntr=P(pp).cntr+1;
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
        %[cmin,xcr,ycr,cmax,xcm,ycm]=get_polydistance(xp,yp,COAST.x_mc,COAST.y_mc,Lcrit);
        [O.timenum_profile] = addVECTOR(O.timenum_profile, timenum1);
        [O.adt_profile] = addVECTOR(O.adt_profile, TIME.adt);
        [O.it_profile] = addVECTOR(O.it_profile, TIME.it);
        [O.c_profile] = addVECTOR(O.c_profile, cmin(:));
        [O.xc_profile] = addVECTOR(O.xc_profile, xcr(:));
        [O.yc_profile] = addVECTOR(O.yc_profile, ycr(:));
        if DUNE.used
            % interpolate dune data to the projection grid (xp,yp)      
            [dmin,xdr,ydr]=get_polydistance(xp,yp,COAST.xdune,COAST.ydune,Lcrit);
            %[dmin,xdr,ydr,dmax,xdm,ydm]=get_polydistance(xp,yp,COAST.xdune,COAST.ydune,Lcrit);
            [O.d_profile] = addVECTOR(O.d_profile, dmin(:));
            [O.xd_profile] = addVECTOR(O.xd_profile, xdr(:));
            [O.yd_profile] = addVECTOR(O.yd_profile, ydr(:));
        end
    end
    
    %% Store all shoreline data to file
    if storedata
        if ~isoctave, warning off, end
        save(fullfile(pwd,O.outputdir,'output.mat'),'O','P','S');
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
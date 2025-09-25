function [O,P] = save_shorelines_netcdf(O, P, TIME, COAST, WAVE, TRANSP, DUNE, MUD)
% function [O,P]=save_shorelines_netcdf(O, P, TIME, COAST, WAVE, TRANSP, DUNE, MUD)
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
%         .x0_mc           : x-coordinates of the initial coastline of all elements at t0
%         .y0_mc           : y-coordinates of the initial coastline of all elements at t0
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
%     DUNE
%         .used            : switch for using the DUNE option 
%     MUD
%         .used            : switch for using the MUD option 
% OUTPUT:
%     O                    : Data-structure with raw instanteneous data, with corresponding fields to the input variables
%         .it              : timestep index (starts at it=0 at model start)
%         .dt              : timestep [years]
%         .tc              : ratio of adaptive time step (tc=0 means using a fixed timestep)
%         .nt              : number of timesteps
%         .timenum         : current model time [days in datenum format]
%     P in NetCDF          : Data structure with projected on a grid which is for the waves and transports also time-averaged data in-between output timesteps, with the following fields:
%         .it              : timestep index (starts at it=0 at model start)
%         .itout           : index for the next output step of the P-structure
%         .cntr            : reset counter for moving average of current output step of the P-structure
%         .tc              : ratio of adaptive time step (tc=0 means using a fixed timestep)
%         .dt              : timestep [years]
%         .nt              : number of timesteps
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

   % Optimized NetCDF writer with unlimited 'space' dimension

   %% Store data at predefined profiles (mainly for plotting)
   [O]=save_profiledata(O,TIME,COAST,DUNE);

   %% Update running-average of physical parameters in P-structure (e.g. average Hs or QS)
   % and store occasionally to the netCDF at each storedata event.
   if isempty(O.storagedate)
      storagedatenum=unique([TIME.timenum0:O.storageinterval:TIME.tend,TIME.tend]);
   else
      storagedatenum=unique(datenum(O.storagedate));
   end
   storedata=0;
   if O.itout<=length(storagedatenum) && TIME.dtsteps==0 || (TIME.tnow+TIME.dt)>=TIME.tend
      storedata=TIME.tnow>=storagedatenum(O.itout) || (TIME.tnow+TIME.dt)>=TIME.tend;
   end

   %% Average data over the current output timestep
   if O.nPspace > 0 && TIME.dtsteps==0
      [O,P] = assemble_projected_data(O, P, COAST, DUNE, MUD, WAVE, TRANSP, TIME,storedata);
   end

   %% === EARLY OUT IF NOT TIME YET ===
   if ~storedata
      return
   end

   %% Open NetCDF
   O.ncid = netcdf.open(O.ncfile, 'NC_WRITE');
   cleanupObj = onCleanup(@() netcdf.close(O.ncid));

   %% === COMPUTE ALONGSHORE DISTANCES AND TIME ===
   TIME.time=TIME.tnow-TIME.timenum0;
   x_mc   = COAST.x_mc(:);
   y_mc   = COAST.y_mc(:);
   dx     = diff(x_mc); dx(isnan(dx)) = 0;
   dy     = diff(y_mc); dy(isnan(dy)) = 0;
   COAST.distx  = [0; cumsum(sqrt(dx.^2 + dy.^2))];
   COAST.distQS = (COAST.distx(1:end-1) + COAST.distx(2:end)) / 2;
   COAST.QS_mc = TRANSP.QS_mc;

   x0_mc   = COAST.x0_mc(:);
   y0_mc   = COAST.y0_mc(:);
   dx1     = diff(x0_mc); dx1(isnan(dx1)) = 0;
   dy1     = diff(y0_mc); dy1(isnan(dy1)) = 0;
   COAST.distx0  = [0; cumsum(sqrt(dx1.^2 + dy1.^2))];
   COAST.distQS0 = (COAST.distx0(1:end-1) + COAST.distx0(2:end)) / 2;
   if DUNE.used==1
   COAST.h0dune_mc = COAST.dcelev_mc - COAST.dfelev_mc;
   end

   %% === ZERO-BASED RECORD INDEX ===
   timeDimID = netcdf.inqDimID(O.ncid,'time');
   [~,  rec]  = netcdf.inqDim(O.ncid,timeDimID);
   spaceDimID        = netcdf.inqDimID(O.ncid,'space');
   [~, existingSpace] = netcdf.inqDim(O.ncid, spaceDimID);

   %% === 1) WRITE TIME SCALARS (always length=1) ===
   for i = 1:size(O.Ivars,1)
       % switch for using/exporting this variable + check if for DUNE or MUD variables
       if O.Ivars{i,1}==1 && (isempty(findstr(O.Ivars{i,2},'dune')) || DUNE.used==1) && (isempty(findstr(O.Ivars{i,2},'mud')) || MUD.used==1)      
           DIM=O.Ivars{i,5};
           STRUCname=O.Ivars{i,6};
           FIELDname=O.Ivars{i,7};
           if strcmpi(STRUCname,'TIME')
               datavals=double(TIME.(FIELDname));
           elseif strcmpi(STRUCname,'COAST')
               datavals=double(COAST.(FIELDname));
           elseif strcmpi(STRUCname,'WAVE')
               datavals=double(WAVE.(FIELDname));
           elseif strcmpi(STRUCname,'STRUC')
               datavals=double(STRUC.(FIELDname));
           else
               datavals=double(eval([STRUCname,'.',FIELDname]));
           end
    
           if DIM==1
               netcdf.putVar(O.ncid, O.ncIDs(i), rec, 1, datavals);
           else
               Ni = numel(datavals);
               Nout = max(Ni, existingSpace);
    
               % build a padded column of length Nout
               col = nan(Nout,1);
               col(1:Ni) = datavals(:);
               netcdf.putVar(O.ncid, O.ncIDs(i), [0 rec], [Nout 1], col);
           end
       end
   end

   %% === 4) WRITE PROJECTED P?STRUCTURES ===
   if O.nPspace > 0
        
        for pp = 1:O.nPspace

            %% Open NetCDF
            P(pp).ncid = netcdf.open(P(pp).ncfile, 'NC_WRITE');
            cleanupObj = onCleanup(@() netcdf.close(P(pp).ncid));
            [~, timLen] = netcdf.inqDim(P(pp).ncid, netcdf.inqDimID(P(pp).ncid, 'time'));        
            spaceLen = length(P(pp).xg);

            for i = 1:size(O.Pvars,1)
                
                % switch for using/exporting this variable + check if for DUNE or MUD variables
                if O.Pvars{i,1}==1 && (isempty(findstr(O.Pvars{i,2},'dune')) || DUNE.used==1) && (isempty(findstr(O.Pvars{i,2},'mud')) || MUD.used==1)      
                    
                    DIM=O.Pvars{i,5};
                    FIELDname=O.Pvars{i,7};
                    
                    [~,~,dimids,~] = netcdf.inqVar(P(pp).ncid, O.ncPInstIDs(pp,i));
                    nDim = numel(dimids);
                                    
                    if nDim == 2
                        % assume dims are [space, time]
                        start = [0,     timLen];
                        count = [spaceLen, 1];
                    elseif nDim == 1
                        [dimName, ~] = netcdf.inqDim(P(pp).ncid, dimids(1));
                        if strcmp(dimName,'time')
                           start = timLen;
                           count = 1;
                        else
                           start = 0;
                           count = spaceLen;
                        end
                    end
                    
                    netcdf.putVar(P(pp).ncid, O.ncPInstIDs(pp,i), start, count, P(pp).(FIELDname));
                end
            end
        end
    end

   %% === ADVANCE COUNTER & FINAL CLOSE ===
   if (O.itout < numel(storagedatenum))
      O.itout = O.itout + 1;
   end
   %netcdf.close(O.ncid);

   if (TIME.tnow + TIME.dt) >= TIME.tend
      clear save_shorelines_netcdf
   end
end

%% Project output data on a grid to P structure (with fixed grid position & mean of transport and waves)
function [O,P] = assemble_projected_data(O, P, COAST, DUNE, MUD, WAVE, TRANSP, TIME, storedata)
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
        if ~isempty(P(pp).xc)
            xc=P(pp).xc;
            yc=P(pp).yc;
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
                P(pp).(fields1{ff})=zeros(length(P(pp).xg),1);
            end
            for ff=1:length(fields2)
                P(pp).(fields2{ff})=zeros(length(P(pp).xg),1);
            end
        end

        % continuously update the data stored as a moving average in the P structure
        wghtNR=[cntr-1,1]/cntr;
        HStdp0=max(P(pp).HStdp,1e-6);
        for ff=1:length(fields1)
            P(pp).(fields1{ff}) = P(pp).(fields1{ff})*wghtNR(1) + scalars.(fields1{ff})(:)*wghtNR(2);
        end
        wghtHS=[HStdp0.^2,scalars.HStdp(:).^2];
        for ff=1:length(fields2)
            sdr = sind(P(pp).(fields2{ff})).*wghtNR(1).*wghtHS(:,1) + sind(vectors.(fields2{ff})(:))*wghtNR(2).*wghtHS(:,2);
            cdr = cosd(P(pp).(fields2{ff})).*wghtNR(1).*wghtHS(:,1) + cosd(vectors.(fields2{ff})(:))*wghtNR(2).*wghtHS(:,2);
            P(pp).(fields2{ff}) = mod(atan2d(sdr,cdr),360);
        end

        % store time variables
        P(pp).tc=TIME.tc;
        P(pp).it=TIME.it;
        P(pp).dt=TIME.dt;
        P(pp).nt=TIME.nt;
        if ind>1 && TIME.tc~=0
            P(pp).dt=(TIME.tnow-P(pp).timenum)/P(pp).cntr/365;
        end
        P(pp).time=TIME.tnow-TIME.timenum0;
        P(pp).timenum=TIME.tnow;

        % store the cross-shore location of the coastline w.r.t. reference
        % line
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
            P(pp).zc=cmin(:);
            P(pp).xc=xcr(:);
            P(pp).yc=ycr(:);
            P(pp).volc=cmin(:).*P(pp).h0;

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
                P(pp).zd=dmin(:);
                P(pp).xd=xdr(:);
                P(pp).yd=ydr(:);
                P(pp).hd0=P(pp).dcelev-P(pp).dfelev;
                P(pp).vold=dmin(:).*P(pp).hd0;
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
function [O]=save_profiledata(O,TIME,COAST,DUNE)
    % these data can be plotted over time
    if ~isempty(O.xyprofiles)
        xp=O.xyprofiles(:,1);
        yp=O.xyprofiles(:,2);
        
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
            %[dmin,xdr,ydr,dmax,xdm,ydm]=get_polydistance(xp,yp,COAST.xdune,COAST.ydune,Lcrit);
            [O.d_profile] = addVECTOR(O.d_profile, dmin(:));
            [O.xd_profile] = addVECTOR(O.xd_profile, xdr(:));
            [O.yd_profile] = addVECTOR(O.yd_profile, ydr(:));
        end
    end
end

%% function to extend the size of a matrix
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
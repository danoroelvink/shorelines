function [COAST,GROYNE,STRUC,TRANSP]=coastline_change(COAST,WAVE,TRANSP,DUNE,MUD,STRUC,GROYNE,TIME,NOUR,FNOUR,CC)
% function [COAST,GROYNE,STRUC,TRANSP]=coastline_change(COAST,WAVE,TRANSP,DUNE,MUD,STRUC,GROYNE,TIME,NOUR,FNOUR,CC)
%
% Updates the coastline and dune position (berm width) at each time step 
% based on transport gradients and nourishment volumes. Both the sand
% and mud transport are accounted for, as well as cross-shore interaction
% with the dunes. Furthermore, the bypassing of sediment at groynes is 
% computed and its subsequent changes in coastline position.
% 
% INPUT:
%   COAST
%       .i_mc           : number of active coastal element 
%       .x              : x-coordinates of cells for considered coastal element
%       .y              : y-coordinates of cells for considered coastal element
%       .n              : number of coastline points of considered coastal element
%       .nq             : number of qs-points of considered coastal element
%       .h0             : active profile height of considered element
%       .ds0            : grid cell size [m]
%       .cyclic         : index whether considered coastal element is cyclical
%       .tanbeta        : slope of the coast [1:n]
%       .PHIc0bnd       : orientation of the shorenormal of the coastline at the boundary points at t0 [°N]
%   <dunes>
%       .qs             : dune erosion volume change that increases the beach width [m3/m/yr]
%       .qss            : dune erosion volume change that increases the beach width, part that is sandy [m3/m/yr]
%       .ql             : dune erosion volume change that does not increase the beach width [m3/m/yr]
%       .qw             : wind transport from the beach to the dune [m3/m/yr]
%       .wberm          : width of the berm/beach [m]
%       .dfelev         : Dune foot elevation [m]
%       .dcelev         : Dune crest elevation [m]
%   <mud>
%       .dndt_mud       : change in mud flat position in [m3/yr/m]
%       .Bf             : mud flat width [m]
%       .Bm             : mangrove width [m]
%       .Bfm            : colonizing mangrove width [m]
%       .dBfdt          : change in mudflat width [m/yr]
%       .dBmdt          : change in mangove width [m/yr]
%       .dBfmdt         : change in colonizing mangrove width [m/yr] 
%   WAVE
%       .diff           : indices of the points with diffracted waves, used for bypassing sediment to.
%   TRANSP
%       .QS             : alongshore transport rates [m3/yr]
%   DUNE
%       .used           : switch for using dunes (0/1)
%   MUD
%       .used           : switch for using mud transport (0/1)
%       .Bmmin          : minimum for the width of the muddy transport zone [m]
%       .Bmmax          : maximum for the width of the muddy transport zone [m]
%   STRUC
%       .shard          : distance along structures [m], with NaNs separating the structures 
%       .xhard          : x-coordinates of the structures [m] (updated)
%       .yhard          : y-coordinates of the structures [m] (updated)
%   GROYNE
%       .n              : number of active groynes
%       .x              : x-coordinates of coastline position at both sides of each structure [Mstruc x 2] [m]
%       .y              : y-coordinates of coastline position at both sides of each structure [Mstruc x 2] [m]
%       .s              : distance along the groyne perimeter for the crossing with the coastline ath both sides of the structure [Mstruc x 2] [m]
%       .QS             : initialization of sediment transport bypass at both sides of the groyne for all groynes [Mx2] at zero [m3/yr]
%       .tipindx        : initialization of the tipindex at both sides of the groyne for all groynes [Mx2]
%       .PHIdn          : difference in angle between groyne and coast at both sides of the groyne for all groynes [Mx2] [°]
%       .strucnum       : index of each groyne inside xhard-yhard, which also contains other structures [M]
%       .idcoast        : index of elements at both sides of the groyne, used for recombining elements after the coastline change [Mx2]
%       .bypassdistpwr  : factor for the redistribution of bypassed sediment over the shadow zone of a groyne (default = 1)
%   TIME
%       .dt             : time step [year]
%   NOUR
%       .growth	        : factor for scaling the nourishment efficiency (from 0 to 1, default is 1)
%       .rate_density   : actual nourishment rate in [m3/m/yr]
%   FNOUR
%       .q_tot          : rate of sediment supply of the shoreface nourishment to the coastline [m3/m/yr]
%   CC
%       .SLRo           : rate of sea level rise [m/yr]
%
% OUTPUT:
%   COAST with updated 'mc' (including MUD and DUNE properties)
%       .dSds_mc        : rate of volume change of the coastline [m3/m/yr]
%       .x_mc           : x-coordinates of coastline points for all coastal elements [m]
%       .y_mc           : y-coordinates of coastline points for all coastal elements [m]
%       .xq_mc          : x-coordinates of qs-points for all coastal elements [m]
%       .yq_mc          : y-coordinates of qs-points for all coastal elements [m]
%       <dunes>
%       .wberm_mc       : width of the berm/beach for all coastal elements [m]
%       .dfelev_mc      : Dune foot elevation for all coastal elements [m]
%       .dcelev_mc      : Dune crest elevation for all coastal elements [m]
%       <mud>
%       .Bf_mc          : mud flat width for all coastal elements [m]
%       .Bm_mc          : mangrove width for all coastal elements [m]
%       .Bfm_mc         : colonizing mangrove width for all coastal elements [m]
%  GROYNE
%       .x              : x-coordinates of coastline position at both sides of each structure [Mstruc x 2] [m]
%       .y              : y-coordinates of coastline position at both sides of each structure [Mstruc x 2] [m]
%       .s              : distance along the groyne perimeter for the crossing with the coastline ath both sides of the structure [Mstruc x 2] [m]
%   STRUC
%       .shard_mc       : distance along structures [m], with NaNs separating the structures 
%       .xhard_mc       : x-coordinates of the structures [m] (updated)
%       .yhard_mc       : y-coordinates of the structures [m] (updated)
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

    eps=1e-3;
    epsgroyne=1;
    idgroyne=[];
    
    %% spatial steps
    n=COAST.n;
    nq=COAST.nq;
    nend=n;
    
    %% Spatial discretisation
    ndev=zeros(1,n);
    COAST.dSds=zeros(1,n);
    COAST.dndt=zeros(1,n);
    COAST.dx=zeros(1,n);
    COAST.dy=zeros(1,n);
       
    %% Determine dominant bypass transport in case of groynes
    if COAST.i_mc==1
        for ig=1:GROYNE.n
            if GROYNE.QS(ig,1)<=0 & GROYNE.QS(ig,2)>=0
            % divergent transports -> set bypass at 0
            GROYNE.QS(ig,:)=0;
            elseif (GROYNE.QS(ig,1)>=0 & GROYNE.QS(ig,2)>=0) || (GROYNE.QS(ig,1)<=0 & GROYNE.QS(ig,2)<=0)    
            % positive or negative transport at both sides -> pick the maximum
            [QSm,imx]=max(abs(GROYNE.QS(ig,:)));
            GROYNE.QS(ig,:)=GROYNE.QS(ig,imx);
            else % GROYNE.QS(ig,1)>=0 & GROYNE.QS(ig,2)<=0
            % converging transport
            GROYNE.QS(ig,:)=sum(GROYNE.QS(ig,:),2);
            end
        end
    end
    
    %% Coastline treatment if a groyne is located at the beginning or end of a coastline section
    xbnd = [nan,nan];
    ybnd = [nan,nan];
    for segmentside=1:2
       if segmentside==2  % in case a groyne is located at the start of a coastal segment
          p1=1;
          p2=2;
          pq=1;
          pq1=1;
          pq2=2;
       elseif segmentside==1  % in case a groyne is located at the end of a coastal segment
          p1=n;
          p2=n-1;
          pq=nq;
          pq1=nq-1;
          pq2=nq;
       end
       if ~isempty(GROYNE.x)
          [idgroyne]=find(GROYNE.idcoast(:,segmentside)==COAST.i_mc);
       end
   
       %% Force the bypass transport rate in the area that is shadowed by the groyne
       % enforce it to use at least the first cell if no shadow zone is present
       % distribute the bypass of the groyne over the shadow zone using a triangular distribution
       % a power can be used to scale the distribution (with more deposition close to the breakwater for pwr>1)
       pwr=GROYNE.bypassdistpwr;
       if ~isempty(idgroyne)
          idbypass=WAVE.diff;
          if ~(sum(idbypass)>=1)
             if segmentside==1
                idbypass(nq)=1;
             else
                idbypass(1)=1;
             end
          end
          if GROYNE.QS(idgroyne,segmentside)<0 && segmentside==1
             TRANSP.QS(nq)=0;
             QSmin=min(TRANSP.QS(idbypass));
             ishadow=find(TRANSP.QS<=max(GROYNE.QS(idgroyne,segmentside),QSmin),1,'last');
             if isempty(ishadow)
                ishadow=nq;
             end
             fraction=max(([1:nq]-ishadow+1)/(nq-ishadow+1),0).^pwr;
             if ishadow<=length(TRANSP.QS)
                TRANSP.QS(ishadow:nq)=TRANSP.QS(ishadow:end)+fraction(ishadow:end).*GROYNE.QS(idgroyne,segmentside);
             end
          elseif GROYNE.QS(idgroyne,segmentside)>0 && segmentside==2
             TRANSP.QS(1)=0;
             QSmax=max(TRANSP.QS(idbypass));
             ishadow=find(TRANSP.QS>=min(GROYNE.QS(idgroyne,segmentside),QSmax),1,'first');
             if isempty(ishadow)
                ishadow=1;
             end
             fraction=max(([nq:-1:1]-(nq-ishadow))/(ishadow),0).^pwr;
             if ishadow>0
                 TRANSP.QS(1:ishadow)=TRANSP.QS(1:ishadow)+fraction(1:ishadow).*GROYNE.QS(idgroyne,segmentside);
             end
          else
              TRANSP.QS(pq)=GROYNE.QS(idgroyne,segmentside);
          end
       end
       
        %% Determine the location where the coastline attaches to the perimeter of a groyne
        if ~isempty(idgroyne)
            % determine the which groyne is used (index 'ist'), the groyne polygon (xhard, yhard) and path along groyne perimeter (shard)
            ist=GROYNE.strucnum(idgroyne);
            [~    ,xhard,~   ,~ , ~ ] = get_one_polygon( STRUC.shard,STRUC.xhard,ist );
            [shard,yhard,~   ,~ , ~ ] = get_one_polygon( STRUC.shard,STRUC.yhard,ist );
            try
                COAST.ds(p1)=hypot(COAST.x(p2)-COAST.x(p1),COAST.y(p2)-COAST.y(p1));
            catch
                ds0=COAST.ds0;
                if ~isscalar(ds0)
                    ds0=min(ds0(:,3));
                end
                COAST.ds(min(max(p1,1),length(COAST.x)))=ds0;
                %figure;plot(COAST.x_mc,COAST.y_mc,'k.-');hold on;plot(STRUC.xhard,STRUC.yhard,'r-');
                fprintf('Please make sure that the groynes are not located exactly at the edge of the model or too close to each other!!! (with just 1 cell elements left inbetween). Model uses reference ds0.\n');
            end
            
            % correct for relative angle of coastline to the groyne. 
            dnratio=1.5./(1-cosd(min(GROYNE.PHIdn(idgroyne,segmentside),180)));
            dnratio=1; %min(dnratio,1);
            
            % compute the coastal change (dn) at the coastline point adjacent to the groyne
            %segmentsign=(segmentside*2-3);
            COAST.dSds(p1)=(TRANSP.QS(pq2)-TRANSP.QS(pq1))/(COAST.ds(p1)/2); 
            COAST.dndt(p1)=-COAST.dSds(p1)/COAST.h0(p1);                          % dndt in [m/yr] = alongshore gradient in sediment transport in [m3/yr/m] divided by 'active height' [m]
            if MUD.used
                COAST.dndt(p1)=COAST.dndt(p1)+COAST.dndt_mud(p1);
            end
            dn=dnratio*(COAST.dndt(p1)+NOUR.growth*NOUR.rate_density(p1)/COAST.h0(p1)-CC.SLRo/COAST.tanbeta+FNOUR.q_tot(p1)/COAST.h0(p1))*TIME.dt;      % dndt in [m/yr] + nourishment in [m3/m/yr] divided by the active height [m] 
            
            % Make sure that the intersection of the coastline with the groyne is not beyond the tip of the groyne
            % GROYNE.s is the distance along the groyne polygon where the crossing is found
            shard_tip=shard(GROYNE.tipindx(idgroyne,segmentside));
            
            if segmentside==2
                if max(GROYNE.s(idgroyne,segmentside)-dn,shard_tip)>max(shard)
                    % extend structure if needed -> so the leg will be extended in landward direction
                    sxt=max(GROYNE.s(idgroyne,segmentside)-dn,shard_tip)-max(shard)+epsgroyne;
                    dx=diff(xhard(end-1:end));
                    dy=diff(yhard(end-1:end));
                    dxt=sxt*dx/((dx.^2+dy.^2).^0.5);
                    dyt=sxt*dy/((dx.^2+dy.^2).^0.5);
                    shard(end)=shard(end)+sxt;
                    xhard(end)=xhard(end)+dxt;
                    yhard(end)=yhard(end)+dyt;
                    [STRUC.xhard,STRUC.yhard]=insert_section(xhard,yhard,STRUC.xhard,STRUC.yhard,ist);
                    [STRUC.shard]=insert_section(shard,STRUC.shard,ist);
                end
                GROYNE.s(idgroyne,segmentside)=min(max(GROYNE.s(idgroyne,segmentside)-dn,shard_tip),max(shard)-eps);
            end
            if segmentside==1
                if min(GROYNE.s(idgroyne,segmentside)+dn,shard_tip)<min(shard)
                    % extend structure if needed -> so the leg will be extended in landward direction
                    sxt=min(GROYNE.s(idgroyne,segmentside)+dn,shard_tip)-min(shard)-epsgroyne;
                    dx=diff(xhard(end-1:end));
                    dy=diff(yhard(end-1:end));
                    dxt=-sxt*dx/((dx.^2+dy.^2).^0.5);
                    dyt=-sxt*dy/((dx.^2+dy.^2).^0.5);
                    shard(1)=shard(1)+sxt;
                    xhard(1)=xhard(1)+dxt;
                    yhard(1)=yhard(1)+dyt;
                    [STRUC.xhard,STRUC.yhard]=insert_section(xhard,yhard,STRUC.xhard,STRUC.yhard,ist);
                    [STRUC.shard]=insert_section(shard,STRUC.shard,ist);
                end
                GROYNE.s(idgroyne,segmentside)=max(min(GROYNE.s(idgroyne,segmentside)+dn,shard_tip),min(shard)+eps);
            end
            
            % find the location where the coastline attaches to the perimeter of the groyne polygon
            GROYNE.x(idgroyne,segmentside)=interp1(shard,xhard,GROYNE.s(idgroyne,segmentside));
            GROYNE.y(idgroyne,segmentside)=interp1(shard,yhard,GROYNE.s(idgroyne,segmentside));
            if ~isnan(GROYNE.x(idgroyne,segmentside))
                xbnd(segmentside)=GROYNE.x(idgroyne,segmentside);
                ybnd(segmentside)=GROYNE.y(idgroyne,segmentside);
            end
        end
    end

    %% COMPUTE RATE OF COASTLINE CHANGE, NOT INCLUDING GROYNE EFFECTS
    % dSds : rate of volumetric change in [m3/yr per meter coastline length]
    % dndt : rate of coastline change in [m/yr]
    for i=1:n
        if COAST.cyclic
            im1=get_mod(i-1,n);
            ip1=get_mod(i+1,n);
            im1q=get_mod(i,nq);
            ip1q=get_mod(i+1,nq);
        else
            im1=max(i-1,1);
            ip1=min(i+1,n);
            im1q=max(i,1);
            ip1q=min(i+1,nq);
        end
        
        COAST.ds(i)=hypot(COAST.x(ip1)-COAST.x(im1),COAST.y(ip1)-COAST.y(im1))/2;
        COAST.ds(i)=max(COAST.ds(i),eps);
        COAST.dSds(i)=(TRANSP.QS(ip1q)-TRANSP.QS(im1q))/COAST.ds(i);
        COAST.dndt(i)=-COAST.dSds(i)/COAST.h0(i);
    end
    COAST.dndtnour=NOUR.rate_density./COAST.h0+FNOUR.q_tot./COAST.h0;

    %% Compute Dune position change
    if DUNE.used
        % Add coastline change due to dune exchange
        % COAST.qss is the sand part of the erosion of the dune
        COAST.dndt=COAST.dndt+(COAST.qss-COAST.qw)./COAST.h0;
        % change berm width when a nourishment is introduced, such that dune foot remains in place
        COAST.wberm=COAST.wberm+COAST.dndtnour*TIME.dt;
        % Make sure wberm has a minimum width
        id1=logical(~TRANSP.shadowS_hD);
        id2=logical(TRANSP.shadowS_hD);
        COAST.wberm(id1)=max(COAST.wberm(id1)+(COAST.dndt(id1)+(COAST.qs(id1)+COAST.ql(id1)-COAST.qw(id1))./(COAST.dcelev(id1)-COAST.dfelev(id1)))*TIME.dt,1);
        COAST.wberm(id2)=max(COAST.wberm(id2)+(COAST.dndt(id2)+(COAST.qs(id2)+COAST.ql(id2)-COAST.qw(id2))./(COAST.dcelev(id2)-COAST.dfelev(id2)))*TIME.dt,1);    
        % Update wberm
        [COAST.wberm_mc,~,~]=insert_props(COAST.wberm,COAST.wberm_mc,COAST.i_mc);
        [COAST.dfelev_mc,~,~]=insert_props(COAST.dfelev,COAST.dfelev_mc,COAST.i_mc);
        [COAST.dcelev_mc,~,~]=insert_props(COAST.dcelev,COAST.dcelev_mc,COAST.i_mc);
    end

    %% Compute Mud coast position change
    if MUD.used
        % Add mud coastline change
        COAST.dndt=COAST.dndt_mud;
        % Update Bf and Bm
        COAST.Bf=COAST.Bf+COAST.dBfdt*TIME.dt;
        COAST.Bfm=COAST.Bfm+(COAST.dBfmdt-COAST.dBmdt)*TIME.dt;
        COAST.Bm=COAST.Bm+COAST.dBmdt*TIME.dt;
        COAST.Bf=max(COAST.Bf,0);
        COAST.Bfm=max(COAST.Bfm,0);
        COAST.Bm=min(COAST.Bm,MUD.Bmmax);
        COAST.Bm=max(COAST.Bm,MUD.Bmmin);
        [COAST.Bf_mc,~,~]=insert_props(COAST.Bf,COAST.Bf_mc,COAST.i_mc);
        [COAST.Bm_mc,~,~]=insert_props(COAST.Bm,COAST.Bm_mc,COAST.i_mc);
        [COAST.Bfm_mc,~,~]=insert_props(COAST.Bfm,COAST.Bfm_mc,COAST.i_mc);
    end
    
    %% Compute Coastline position change (dn)
    for i=1:n
        if COAST.cyclic
            im1=get_mod(i-1,n);
            ip1=get_mod(i+1,n);
        else
            im1=max(i-1,1);
            ip1=min(i+1,n);
        end
        dn=(COAST.dndt(i)+COAST.dndtnour(i)-CC.SLRo/COAST.tanbeta)*TIME.dt;
       
        if COAST.preserveorientation==1
            if length(COAST.PHIcxy0)~=length(COAST.x)
                fprintf('Warning : grid length has changed');
            end
            COAST.dx(i)= dn*sind(COAST.PHIcxy0(i)); % no minus sign here
            COAST.dy(i)= dn*cosd(COAST.PHIcxy0(i)); 
        elseif i==1
            COAST.dx(i)= dn*sind(COAST.PHIcxy0(1)); % no minus sign here
            COAST.dy(i)= dn*cosd(COAST.PHIcxy0(1));           
        elseif i==n
            COAST.dx(i)= dn*sind(COAST.PHIcxy0(end)); % no minus sign here
            COAST.dy(i)= dn*cosd(COAST.PHIcxy0(end)); 
        else %if COAST.cyclic 
            COAST.dx(i)=-dn*(COAST.y(ip1)-COAST.y(im1))/(2*COAST.ds(i));
            COAST.dy(i)= dn*(COAST.x(ip1)-COAST.x(im1))/(2*COAST.ds(i));
        end
    end
    
    %% Update coastline positions, all points
    x0=COAST.x(:);
    y0=COAST.y(:);
    for i=1:n
        COAST.x(i)=x0(i)+COAST.dx(i);
        COAST.y(i)=y0(i)+COAST.dy(i);
    end
    if COAST.cyclic
        x1=(COAST.x(n)+COAST.x(1))/2;
        y1=(COAST.y(n)+COAST.y(1))/2;
        COAST.x(1)=x1;
        COAST.y(1)=y1;
        COAST.x(n)=x1;
        COAST.y(n)=y1;
    end
    
    %% Overwrite coastline position if groyne at begin of coastline section
    if ~isnan(xbnd(2))
        COAST.x(1)=xbnd(2);
        COAST.y(1)=ybnd(2);
    end
    %% Overwrite coastline position if groyne at end of coastline section
    if ~isnan(xbnd(1))
        COAST.x(n)=xbnd(1);
        COAST.y(n)=ybnd(1);
    end
    
    % Replace nans with -1e10 in this way the number of elements stays n_mc
    if isempty(find(~isnan(COAST.x)))
        COAST.x=-1e10;
        COAST.y=-1e10;
        COAST.dSds=-1e10;
    end

    % Insert dSds in dSds_mc
    if  COAST.i_mc==1
        [COAST.dSds_mc]=COAST.dSds;
    else
        [COAST.dSds_mc]=[COAST.dSds_mc,nan,COAST.dSds];
    end
    
    %% insert new section into x_mc and y_mc
    [COAST.x_mc,COAST.y_mc]=insert_section(COAST.x,COAST.y,COAST.x_mc,COAST.y_mc,COAST.i_mc);
    
    %% insert new section TRANSP.QS into TRANSP.QS_mc
    [TRANSP.QS_mc]=insert_section(TRANSP.QS,TRANSP.QS_mc,COAST.i_mc);

    %% make transport points xq_mc and yq_mc
    [COAST]=get_transportpoints(COAST,COAST.i_mc);
    
    % Check area
    %A1 = polyarea(x0(:),y0(:));
    %A2 = polyarea(COAST.x(:),COAST.y(:));
    %A3 = sum(COAST.ds_mc(~isnan(COAST.ds_mc)).*COAST.y_mc(~isnan(COAST.ds_mc)))
    %fprintf('Area=%12.0f, Area2=%12.0f      ',A1,A2);
end

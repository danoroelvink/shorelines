function [ COAST,SPIT ] = find_overwash_mc(COAST,WAVE,SPIT,STRUC,TRANSP,DUNE,MUD,TIME)
% function [ COAST,SPIT ] = find_overwash_mc(COAST,WAVE,SPIT,STRUC,TRANSP,DUNE,MUD,TIME)
%
% The overwash is computed from the seaward side of barriers to the landward side. 
% Make sure to select coastline points that are not in shadow of coastline or structures.
% The model will for each coastline point construct a transect in the direction of the waves. 
% Intersections of this ray with other coastal elements are checked. If found, then 
% the distance to the element on the backside of the barrier is computed, 
% which may trigger overwash. The sediment is distributed over points at the backside
% of the barrier. 
%
% INPUT: 
%    COAST
%        .x_mc        : x-coordinates of grid [m]
%        .y_mc        : y-coordinates of grid [m]
%        .n_mc        : number of coastal elements
%        .cyclic_mc   : index of cyclicality for each coastal element
%        .shadow_mc   : index of coastline points with a shadow
%        .ql_mc       : dune parameter
%        .dcelev_mc   : dune crest elevation w.r.t. MSL (m) (only used if ldbdune is empty)
%    TIME
%        .dt          : timestep of the model [yrs]
%    SPIT
%        .method      : overwash method ('default' or 'dune_overwash')
%        .owtimescale : timescale for overwash (i.e. what part of the deficit is moved to the backbarrier)
%        .spitwidth   : equilibrium width of the spit [m]
%        .bheight     : berm height used for overwash funciton (i.e. added to Dsf or Dbb)
%        .Dbb         :	underwater part of active height for shoreface [m]
%        .Dsf         : underwater part of active height for back-barrier [m]
% 
% OUTPUT:
%    COAST
%        .x_mc        : x-coordinates of grid [m] (corrected for overwash)
%        .y_mc        : x-coordinates of grid [m] (corrected for overwash)
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

    %% INITIALIZE VARIABLES
    eps=5;
    yesplot=0;
    dxt = zeros(size(COAST.x_mc));
    dyt = zeros(size(COAST.x_mc));
    
    %% For each coast section
    % In the following, ii denotes the section from where we check the spit width
    % and jj the section of the point on the other side of the spit
    for ii=1:COAST.n_mc
        
        [COAST,WAVE,~]=get_segmentdata(COAST,WAVE,TRANSP,DUNE,MUD,ii);
        if length(COAST.x)<3
            return
        end
        SPIT.width=zeros(size(COAST.x));
        
        % compute phiw_coast at nodes from PHIWI at transport points
        phiw_coast=mod(atan2d(sind(0.5*(WAVE.PHI(1:end-1)+WAVE.PHI(2:end))),cosd(0.5*(WAVE.PHI(1:end-1)+WAVE.PHI(2:end)))),360.0);
        
        % loop over coastal elements
        for i=1:length(COAST.x)
            
            if ~COAST.shadow(i)
                if COAST.cyclic
                    im1=i-1;if im1<1;im1=COAST.n-1;end
                    ip1=i+1;if ip1>COAST.n;ip1=2;end
                else
                    im1=max(i-1,1);
                    ip1=min(i+1,COAST.n);
                end
                
                % construct line in coast normal direction x_nl, y_nl (normal landward)
                phic=360-atan2d(COAST.y(ip1)-COAST.y(im1),COAST.x(ip1)-COAST.x(im1));
                x_nl=[COAST.x(i)-eps*sind(phic),COAST.x(i)-SPIT.spitwidth*sind(phic)];
                y_nl=[COAST.y(i)-eps*cosd(phic),COAST.y(i)-SPIT.spitwidth*cosd(phic)];
                
                % get intersections with x_mc, y_mc
                [xcr,ycr,~,~,ind_mc,~,ui,~]=get_intersections(COAST.x_mc,COAST.y_mc,x_nl,y_nl);
                
                if ~isempty(xcr)&~isnan(xcr)&~isnan(ycr)
                    % find closest point, distance, reference points and weights
                    [distance,closest] = min(hypot(xcr-COAST.x(i),ycr-COAST.y(i)));
                    SPIT.width(i)      = distance;
                    j                  = ind_mc(closest);
                    weight_j           = 1-ui(closest);
                    weight_jp1         = ui(closest);
                    if weight_j==1
                        jp1=j;
                        weight_jp1=1;
                    else
                        jp1=j+1;
                    end
                    
                    % compute shift of coastline due to overwashing, depending on method
                    switch lower(SPIT.method)
                        case ('default')
                            if SPIT.owtimescale>0
                                dn_sea  = -TIME.dt/SPIT.owtimescale*(SPIT.spitwidth-SPIT.width(i))*max(cosd(phiw_coast(i)-phic),0);
                                dn_back = -dn_sea  *(SPIT.Dsf+SPIT.bheight)/(SPIT.Dbb+SPIT.bheight);
                            else
                                dn_sea  = -SPIT.owscale*(SPIT.spitwidth-SPIT.width(i))*max(cosd(phiw_coast(i)-phic),0);
                                dn_back = -dn_sea  *(SPIT.Dsf+SPIT.bheight)/(SPIT.Dbb+SPIT.bheight);
                            end
                        case ('dune_overwash')
                            dn_sea  = 0; % Already taken care of in coastline_change
                            dn_back = TIME.dt*COAST.ql(i)/(COAST.dcelev(i)+SPIT.Dbb);                     
                    end
                    
                    
                    % apply shift to seaward point
                    ds_sea=hypot(COAST.x(ip1)-COAST.x(im1),COAST.y(ip1)-COAST.y(im1));
                    dx=-dn_sea*(COAST.y(ip1)-COAST.y(im1))/ds_sea;
                    dy= dn_sea*(COAST.x(ip1)-COAST.x(im1))/ds_sea;
                    %
                    % update coastline change in global arrays
                    dxt(i+COAST.i1-1)=dxt(i+COAST.i1-1)+dx;
                    dyt(i+COAST.i1-1)=dyt(i+COAST.i1-1)+dy;
                    
                    % apply shift to backbarrier point
                    try
                        ds_back=hypot(COAST.x_mc(j+1)-COAST.x_mc(j),COAST.y_mc(j+1)-COAST.y_mc(j));
                    catch
                        disp('problem')
                    end
                    %dn_back=dn_back*ds_sea/2/ds_back;  % ds_back taken only over 1 cell
                    
                    dxt(j)=-dn_back*(COAST.y_mc(j+1)-COAST.y_mc(j))/ds_back*weight_j;
                    dyt(j)= dn_back*(COAST.x_mc(j+1)-COAST.x_mc(j))/ds_back*weight_j;
                    dxt(jp1)=-dn_back*(COAST.y_mc(jp1)-COAST.y_mc(j))/ds_back*weight_jp1;
                    dyt(jp1)= dn_back*(COAST.x_mc(jp1)-COAST.x_mc(j))/ds_back*weight_jp1;
                end
            end
        end
    end
    
    %% Check cyclic boundaries
    for ii=1:COAST.n_mc
        [COAST,WAVE,TRANSP,i1,i2] = get_segmentdata(COAST,WAVE,TRANSP,DUNE,MUD,ii);
        if COAST.cyclic
            if hypot(dxt(i1)-dxt(i2),dyt(i1)-dyt(i2)) > 1 & COAST.shadow(1)
                % back barrier side - could have contributions from both sides of cyclic point
                sumdxt=dxt(i1)+dxt(i2);
                sumdyt=dyt(i1)+dyt(i2);
                dxt(i1)=sumdxt;
                dxt(i2)=sumdxt;
                dyt(i1)=sumdyt;
                dyt(i2)=sumdyt;
            end
        end
    end
    
    %% Now update the coastlines
    COAST.x_mc = COAST.x_mc + dxt;
    COAST.y_mc = COAST.y_mc + dyt;
    
    
    if yesplot
        plot(COAST.x_mc,COAST.y_mc,'--k')
        hold off
        drawnow
    end
end

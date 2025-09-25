function [COAST,STRUC,GROYNE]=prepare_grid_groyne(COAST,STRUC,yesplot)
% function [COAST,STRUC,GROYNE]=prepare_grid_groyne(COAST,STRUC,yesplot)
% 
% The routine initializes the coastline grid, taking into account groynes 
% intersecting with the coastline. 
% - Offshore breakwaters, revetments and groynes are differentiated (idgroyne)
%   with an index 0 for revetments and offshore breakwaters and for  
%   groynes an index number (larger than 0)
% - The groyne order of coordinates is flipped if it is not clockwise. 
%   The xhard and yhard are then updated. 
% - The distance along the perimeter
%   of the structures is computed (shard). 
% - The model will find intersections between groynes and the coast. 
%   If necessary the ends of the groyne are extended a bit virtually to
%   allow for structures at the boudnary. 
% - if a cyclical element (e.g. an offshore breakwater) intersects the coastline
%   then it is checked if the first and last point can be moved 'inland', 
%   basically re-organizing the structure (xhard and yhard) order of coordinates
%   such that it can act as a groyne afterwards (e.g. when a tombolo has developed). 
% - For the groynes their orientation (phigroyne), orientation of groyne at 
%   location where it crosses the coast (phigroynesize) and the local coastline 
%   orientation at these crossing locations are determined (phicoast), 
%   as well as the difference in orientation with respect to the coast (PHIdn), 
%   which is relevant for the wave diffraction later on. 
% - The coastline will be split at the groynes (x_mc and y_mc), and the 
%   index of the coastline sections at both sides of each groyne stored (idcoast)
%   as well as the index of groynes at the ends of the coastline sections (i_groyne)
%   These indices are relevant later on when the coastal elements are recombined.
%
% INPUT:  
%    COAST
%        .x_mc                    : x-coordinates of the coastline [m]
%        .y_mc                    : y-coordinates of the coastline [m]
%        .ds0                     : grid cell size [m]
%        .PHIcxy_mc               : (optional) orientation at coastline points [°N]
%    STRUC
%        .bypassdistpwr           : factor for the redistribution of bypassed sediment over the shadow zone of a groyne (default = 1)
%        .xhard                   : x-coordinates of the structures [m]
%        .yhard                   : y-coordinates of the structures [m]
%        .nrevet                  : number of revetments (i.e. the last structures  in xhard/yhard are these revetments)
%    yesplot                      : switch for debug plotting (0/1)
% 
% OUTPUT:
%    STRUC
%        .nhard                   : number of structures (updated)
%        .xhard                   : x-coordinates of the structures [m] (updated)
%        .yhard                   : y-coordinates of the structures [m] (updated)
%        .idgroyne                : index for each structure, indicating whether it is a groyne (with value = groyne number) or offshore breakwater/revetment (with value = 0) . 
%        .shard                   : distance along structures [m], with NaNs separating the structures 
%        .PHIs                    : orientation of path along the structure (xhard and yhard) [°N], with NaNs separating the structures
%    GROYNE
%        .n                       : number of active groynes
%        .x                       : x-coordinates of crossings of each structure with the coast [Mx2] with at each line the first and second crossing [m]
%        .y                       : y-coordinates of crossings of each structure with the coast [Mx2] with at each line the first and second crossing [m]
%        .s                       : distance along the groyne perimeter with the crossing with the coast [Mx2] with at each line the first and second crossing [m]
%        .xg                      : x-coordinate of the point inside the groyne (one value for each groyne) [m]
%        .yg                      : y-coordinate of the point inside the groyne (one value for each groyne) [m]
%        .QS                      : initialization of sediment transport bypass at both sides of the groyne for all groynes [Mx2] at zero [m3/yr]
%        .tipindx                 : initialization of the tipindex at both sides of the groyne for all groynes [Mx2]
%        .phicoast                : coastline orientation at crossing point of groyne and coast for all groynes [Mx2] [°N]
%        .phigroyne               : average orientation of the groyne (i.e. begin and end point compared to average of rest of groyne) for all groynes [M] [°N]
%        .phigroyneside           : orientation of the groyne itself at both crossings with the coast for all groynes [Mx2] [°N]
%        .phicoastside            : coastline orientation at left and right side of groyne for all groynes [Mx2] [°N]
%        .PHIdn                   : difference in angle between groyne and coast at both sides of the groyne for all groynes [Mx2] [°]
%        .strucnum                : index of each groyne inside xhard-yhard, which also contains other structures [M]
%        .idcoast                 : index of segments at both sides of the groyne, used for recombining elements after the coastline change [Mx2]
%        .bypassdistpwr           : factor for the redistribution of bypassed sediment over the shadow zone of a groyne (default = 1)
%    COAST
%        .x_mc                    : x-coordinates of the coastline elements after splitting at groynes [m]
%        .y_mc                    : y-coordinates of the coastline elements after splitting at groynes [m]
%        .n_mc                    : number of coastline elements after splitting at groynes (updated)
%        .x_mc0gridgroyne         : x-coordinates of the coastline before 'prepare groyne' [m]
%        .y_mc0gridgroyne         : y-coordinates of the coastline before 'prepare groyne' [m]
%        .i_groyne                : stores the index of the groynes (or locations of the splitting) along the coastline grid x_mc/y_mc (e.g. three groynes and four coastal sections : [0 0 0 0 0 0 1 0 0 0 0 0 2 0 0 0 3 0 0])
%        .BNDgroyne               : index showing the presence of structures at the start-end of each coastal section, and their structure id [number of coastal elements x 2-sides]
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
    
    %fprintf('  Initialize groynes (grid) \n');   
    %% insert  points at the intersection between structure and coastline if exists  ( Added by M.Ghonim )
    GROYNE.s=[];
    GROYNE.x=[];
    GROYNE.y=[];
    GROYNE.QS=[];
    GROYNE.n=0;
    GROYNE.strucnum=[];
    GROYNE.idcoast=[];
    GROYNE.bypassdistpwr=STRUC.bypassdistpwr;
    STRUC.idgroyne=[];
    COAST.BNDgroyne=zeros(COAST.n_mc,2);
    COAST.x_mc0gridgroyne=COAST.x_mc;
    COAST.y_mc0gridgroyne=COAST.y_mc;
    GROYNE.xg=[];
    GROYNE.yg=[];
    
    %% Find intersections with groynes, split coastline and store groyne intersections
    if length(STRUC.xhard)>1
        n_st=length(find(isnan(STRUC.xhard)))+1;                                    % Number of groynes
        STRUC.nhard=n_st;
        ii=0;
        for ist=1:n_st
            [ xs,ys,~,~,~ ] = get_one_polygon( STRUC.xhard,STRUC.yhard,ist);xs=xs(:)';ys=ys(:)';
            STRUC.idgroyne(ist)=0;

            % compute the distance (ss) along the polygon of the hard structure, and add it to 'STRUC.shard'.
            ss=[0,cumsum((diff(xs).^2+diff(ys).^2).^0.5)];
            if ist==1
                STRUC.shard=[ss];
            else
                STRUC.shard=[STRUC.shard,nan,ss];
            end
            % orientation of groyne structure lines
            STRUC.PHIs=mod(atan2d(diff(STRUC.xhard),diff(STRUC.yhard)),360);
            
            % a groyne strcuture should at least be 2 points long, and it should not be a revetment
            if length(xs)>2 && ist<=n_st-STRUC.nrevet 
                 % find intersections of current structure with full coastline
                [str_int_x,str_int_y,indc,inds]=get_intersections(COAST.x_mc,COAST.y_mc,xs,ys);
                
                % make cyclical structures clockwise
                scyclic=get_cyclic(xs,ys,COAST.ds0);
                if scyclic==1
                    cw=get_clockpoly(xs,ys);
                    if cw~=1
                        xs=fliplr(xs);
                        ys=fliplr(ys);
                        % add mirrored xs and ys to the xhard and yhard
                        [STRUC.xhard,STRUC.yhard]=insert_section(xs,ys,STRUC.xhard,STRUC.yhard,ist);
                        [str_int_x,str_int_y,indc,inds]=get_intersections(COAST.x_mc,COAST.y_mc,xs,ys);
                        % compute the distance (ss) along the polygon of the hard structure, and add it to 'STRUC.shard'.
                        ss=fliplr(ss(end)-ss);
                        % add new distance along structure to shard
                        [STRUC.shard]=insert_section(ss,STRUC.shard,ist);
                        % orientation of groyne structure lines
                        STRUC.PHIs=mod(atan2d(diff(STRUC.xhard),diff(STRUC.yhard)),360);
                    end
                end
                
                % groyne defined at boundary then check if it crosses at x(1)-ds or at x(end)+ds
                if length(str_int_x)==1 
                    % check if not cyclic
                    dist=((COAST.x_mc(end)-COAST.x_mc(1)).^2+(COAST.y_mc(end)-COAST.y_mc(1)).^2).^0.5;
                    if dist>COAST.ds0
                        % compute distance of groyne to boundary on left or right side
                        dist1=hypot(COAST.x_mc(1)-str_int_x,COAST.y_mc(1)-str_int_y);
                        dist2=hypot(COAST.x_mc(end)-str_int_x,COAST.y_mc(end)-str_int_y);
                        
                        % extend grid at x1 or xend
                        dx=diff(COAST.x_mc);
                        dy=diff(COAST.y_mc);
                        L=(dx.^2+dy.^2).^0.5;
                        dxs=10*COAST.ds0.*dx./L;
                        dys=10*COAST.ds0.*dy./L;
                        x_mcext=COAST.x_mc;
                        y_mcext=COAST.y_mc;
                        if dist1<2*COAST.ds0 
                            x_mcext=[COAST.x_mc(1)-dxs(1),COAST.x_mc];
                            y_mcext=[COAST.y_mc(1)-dys(1),COAST.y_mc];
                        elseif dist2<2*COAST.ds0
                            x_mcext=[COAST.x_mc,COAST.x_mc(end)+dxs(end)];
                            y_mcext=[COAST.y_mc,COAST.y_mc(end)+dys(end)];
                        end 
                        
                        % find intersections with groynes for the extended grid
                        [str_int_x2,str_int_y2,indc,inds]=get_intersections(x_mcext,y_mcext,xs,ys);
                        if length(str_int_x2)>1
                            % reconstruct grid which is extended to point of cross-section with breakwater
                            if dist1<COAST.ds0 
                                x_mcbnd=[str_int_x2(find(indc<2,1,'last')),COAST.x_mc];
                                y_mcbnd=[str_int_y2(find(indc<2,1,'last')),COAST.y_mc];
                                indc=indc-1;
                            elseif dist2<COAST.ds0 
                                x_mcbnd=[COAST.x_mc,str_int_x2(find(indc,1,'last'))];
                                y_mcbnd=[COAST.y_mc,str_int_y2(find(indc,1,'last'))];
                                indc(find(indc==max(indc)))=1e10;
                            end
                            str_int_x=str_int_x2;
                            str_int_y=str_int_y2;
                        end
                    end
                end
                
                % flip groyne if needed (i.e. when the first crossing is on a higher coastline grid index than the second crossing)
                % or when it is a closed element with the start and end point offshore instead of on land. 
                [~,idc]=sort(indc);
                if ~isempty(idc)
                    if idc(1)>idc(end)
                        if scyclic==1
                            % change start-end location, such that is located inland
                            idcst = find([1:length(ss)]>inds(1) & [1:length(ss)]<=inds(2));
                            if ~isempty(idcst)
                                % in case there are a few points inside the coast
                                idcst=idcst(round(mean([1:length(idcst)])));
                                xs=[xs(idcst:end-1),xs(1:idcst)];
                                ys=[ys(idcst:end-1),ys(1:idcst)];
                            else
                                % in case there is not yet a point of the hard structure beyond the coastline
                                idcst1 = find([1:length(ss)]<inds(1));
                                idcst2 = find([1:length(ss)]>=inds(2));
                                x0=mean(str_int_x); %(xs(idcst1(end))+xs(idcst2(1)))/2;
                                y0=mean(str_int_y); %(ys(idcst1(end))+ys(idcst2(1)))/2;
                                xs=[x0,xs(idcst2(1:end-1)),xs(idcst1),x0];
                                ys=[y0,ys(idcst2(1:end-1)),ys(idcst1),y0];
                            end
                        else
                            % flip groyne if it is an open element that is placed counter-clockwise
                            xs=fliplr(xs);
                            ys=fliplr(ys);
                        end
                        % add mirrored xs and ys to the xhard and yhard
                        [STRUC.xhard,STRUC.yhard]=insert_section(xs,ys,STRUC.xhard,STRUC.yhard,ist);
                        [str_int_x,str_int_y,indc,inds]=get_intersections(COAST.x_mc,COAST.y_mc,xs,ys);
                        % compute the distance (ss) along the polygon of the hard structure, and add it to 'STRUC.shard'.
                        ss=[0,cumsum((diff(xs).^2+diff(ys).^2).^0.5)];
                        % add new distance along structure to shard
                        [STRUC.shard]=insert_section(ss,STRUC.shard,ist);
                        % orientation of groyne structure lines
                        STRUC.PHIs=mod(atan2d(diff(STRUC.xhard),diff(STRUC.yhard)),360);
                    end
                end
                
                % correct the crossing points of the groyne with the coastline if it detects some twice (then use first and last one)
                [~,idu]=unique(str_int_x+i*str_int_y);idu=idu(:)';
                if length(str_int_x)>2
                    idu1=find(indc==min(indc));
                    idu2=find(indc==max(indc));
                    idu=[idu1,idu2]; % also make sure no more than two crossingpoints are used
                    [~,idc2]=sort(indc(idu));
                    idu=idu(idc2);
                    indc=indc(idu);
                    inds=inds(idu);
                    str_int_x=str_int_x(idu);
                    str_int_y=str_int_y(idu);
                elseif length(idu)==1 && length(indc)>=2 
                    % situation where the groyne polygon crosses twice, but exactly at the same grid point. 
                    % Then use the next coastal grid point for the second crossing.
                    % this is not yet used! (since length(idu) should be 2 later on)
                    indc=int32([indc(1)+[0,1]]);
                    if length(unique(inds))==1
                        inds=int32([inds(1)+[0,1]]);
                    end
                    str_int_x(2)=COAST.x_mc(indc(2));
                    str_int_y(2)=COAST.y_mc(indc(2));
                end
                
                % a groyne needs at least 2 unique crossing points with the coastline
                % determine the x and y location where the groyne attaches to the coast (in GROYNE.x and GROYNE.y)
                % and store the coastal segment numbers the groyne is connected to (in GROYNE.idcoast)
                % then cut the coast in two pieces at both sides of the groyne. 
                if length(idu) == 2 
                    % index linking the hard structure # (ist) to the groyne # (ii)
                    ii=ii+1;
                    STRUC.idgroyne(ist)=ii;
                    
                    % valid groyne, with two intersections
                    GROYNE.n=GROYNE.n+1;
                    GROYNE.strucnum(ii)=ist;
                    GROYNE.s(ii,:)=interp1([1:length(ss)],ss,[inds(1),inds(end)]);
                    GROYNE.x(ii,:)=str_int_x+diff(str_int_x)*1e-5*[1,-1];
                    GROYNE.y(ii,:)=str_int_y+diff(str_int_y)*1e-5*[1,-1];
                    indc0=[];
                    if indc(1)~=1e10 && indc(2)~=1e10
                        in=find(inpolygon(COAST.x_mc,COAST.y_mc,xs,ys));
                        if ~isempty(in)
                            indc0=round(mean(in));
                            indc1=round(min(in))-1;
                            indc2=round(max(in))+1;
                        else
                            indc1=max(floor(indc(1)),1);
                            indc2=min(ceil(indc(2)),length(COAST.x_mc));
                        end  
                    elseif isfield(COAST,'PHIcxy_mc')
                        if indc(1)==1e10 
                            indc1=1;
                            indc2=1;
                        end
                        if indc(2)==1e10 
                            indc1=length(COAST.PHIcxy_mc);
                            indc2=indc1;
                        end
                    end
                    if ~isempty(indc0)
                        GROYNE.xg(ii)=COAST.x_mc(indc0);  % xg and yg indicate the point inside the groyne.
                        GROYNE.yg(ii)=COAST.y_mc(indc0);
                    else
                        GROYNE.xg(ii)=mean(GROYNE.x(ii,:));
                        GROYNE.yg(ii)=mean(GROYNE.y(ii,:));
                    end
                    GROYNE.QS=zeros(GROYNE.n,2);
                    GROYNE.tipindx=zeros(GROYNE.n,2);
                    GROYNE.phicoast(ii,1)=mod(360-atan2d(diff(str_int_y),diff(str_int_x)),360);
                    dxgroyne=mean(xs(2:end-1))-(xs(1)+xs(end))/2;
                    dygroyne=mean(ys(2:end-1))-(ys(1)+ys(end))/2;
                    GROYNE.phigroyne(ii,1)=mod(90-atan2d(dygroyne,dxgroyne),360);
                    
                    % compute parameters that are used to determine angle of groyne with respect to the coastline
                    PHIs=mod(atan2d(diff(xs),diff(ys)),360);
                    GROYNE.phigroyneside(ii,1)=mod(PHIs(max(floor(inds(1)),1)),360); %mod(atan2d(str_int_x(1)-GROYNE.xg(ii),str_int_y(1)-GROYNE.yg(ii)),360);
                    GROYNE.phigroyneside(ii,2)=mod(PHIs(min(floor(inds(2)),length(PHIs)))-180,360); %mod(atan2d(str_int_x(2)-GROYNE.xg(ii),str_int_y(2)-GROYNE.yg(ii)),360);
                    if isfield(COAST,'PHIcxy_mc') 
                        GROYNE.phicoastside(ii,1)=COAST.PHIcxy_mc(min(max(indc1,1),length(COAST.PHIcxy_mc))); % coastline orientation at left side of groyne
                        GROYNE.phicoastside(ii,2)=COAST.PHIcxy_mc(min(max(indc2,1),length(COAST.PHIcxy_mc))); % coastline orientation at right side of groyne
                    else
                        GROYNE.phicoastside(ii,1)=GROYNE.phigroyneside(ii,1);
                        GROYNE.phicoastside(ii,2)=GROYNE.phigroyneside(ii,2);
                    end
                    GROYNE.PHIdn(ii,:)=abs(GROYNE.phigroyneside(ii,:)+[90,-90]-GROYNE.phicoastside(ii,:));
                    
                    % find segments index (i_mc) at both sides of the groyne split
                    idnotnans=find(~isnan(COAST.x_mc(1:floor(min(indc)))));  % find the number of segments before        
                    nrsegmentsbefore=length([idnotnans(diff(idnotnans)~=1)])+1; % count the number of segments 
                    
                    % add index of segments at both sides of the groyne
                    GROYNE.idcoast(ii,:)=[nan,nan];
                    
                    % special case where the groyne is at the beginning of the grid
                    % special case where the groyne is at the end of the grid
                    if indc(1)<1
                        GROYNE.idcoast(ii,:)=[nan,1];
                    elseif indc(2)==1e10
                        GROYNE.idcoast(ii,:)=[nrsegmentsbefore,nan];
                    else                      
                        GROYNE.idcoast(ii,:)=nrsegmentsbefore+[0,1];
                    end
                    
                    % make sure the segment-ids of earlier defined groyne-segment splits are adjusted (for segments after the current split)
                    if indc(1)>=1
                    idelementsaftersplit=setdiff(find(GROYNE.idcoast(:,1)>=nrsegmentsbefore),ii);
                    GROYNE.idcoast(idelementsaftersplit,:)=GROYNE.idcoast(idelementsaftersplit,:)+1;  
                    end
                    
                    % add nans to coastline
                    if floor(min(indc))~=0 && ceil(max(indc))~=1e10
                        COAST.x_mc = [COAST.x_mc(1:floor(min(indc))), GROYNE.x(ii,1),NaN,GROYNE.x(ii,2),COAST.x_mc(ceil(max(indc)):end)];
                        COAST.y_mc = [COAST.y_mc(1:floor(min(indc))), GROYNE.y(ii,1),NaN,GROYNE.y(ii,2),COAST.y_mc(ceil(max(indc)):end)];
                    elseif floor(min(indc))==0
                        %COAST.x_mc = [GROYNE.x(ii,2),COAST.x_mc(ceil(max(indc)):end)];
                        %COAST.y_mc = [GROYNE.y(ii,2),COAST.y_mc(ceil(max(indc)):end)];
                    elseif ceil(max(indc))==1e10
                        %COAST.x_mc = [COAST.x_mc(1:floor(min(indc))), GROYNE.x(ii,1)];
                        %COAST.y_mc = [COAST.y_mc(1:floor(min(indc))), GROYNE.y(ii,1)];
                    end
                        
                    % add groyne number to array to make sure that the splitting nans of groynes can be traced back
                    COAST.i_groyne(floor(indc(1))+1)=ist;  
                                                
                elseif length(str_int_x)<2
                    %disp('Less than 2 crossing points of the coast and groyne are found!');
                    %disp('Structure is assumed to be burried below the sand or acts as an offshore breakwater.');
                else
                    disp(num2str(str_int_x))
                end 
            end
        end

        % determine number of coastline segments
        nans=find(isnan(COAST.x_mc));
        COAST.n_mc=length(nans)+1;
        COAST.BNDgroyne=zeros(COAST.n_mc,2);

        if yesplot
           figure(11);
        end

        % Stores the groyne-ids at the position of the open boundaries of the models
        % the indices in COAST.BNDgroyne are the index of the groynes:
        % e.g. [ 0 1; 
        %        1 0; 
        %        0 2; 
        %        2 0]; 
        % This means that there are two groynes. 
        % The first groyne between coastline elements 1 (right side) and 2 (left side). 
        % And the second groyne between element 3 (right side) and element 4 (left side).
        for gg=1:size(GROYNE.idcoast,1)
            if ~isnan(GROYNE.idcoast(gg,1))
            COAST.BNDgroyne(GROYNE.idcoast(gg,1),2)=gg;
            end
            if ~isnan(GROYNE.idcoast(gg,2))
            COAST.BNDgroyne(GROYNE.idcoast(gg,2),1)=gg;
            end
        end
    end
end
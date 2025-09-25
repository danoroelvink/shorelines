function [COAST]=make_sgrid_mc(COAST,TIME,i_mc)
% function [COAST]=make_sgrid_mc(COAST,TIME,i_mc)
%
% The model grid is made in this function. There are two
% methods for creating the grid:
%   griddingmethod 1 : split or merge only the cells that are 
%                      too small or too large. This approach 
%                      reduces diffusion in the grid generation. 
%   griddingmethod 2 : re-interpolates the whole grid if a single
%                      grid cell is too large or too small.
%                      This very stable but does effectively 
%                      creates a small bit of diffusion. 
%
% INPUT: 
%     COAST
%         .x_mc           : x-coordinates of all coastal sections [m]
%         .y_mc           : y-coordinates of all coastal sections [m]
%         .ds0            : preferred grid cell size [m]
%         .griddingmethod : method of grid generation, either:
%                           1=locally split/merge cells if needed, or 
%                           2=re-generate a full smooth grid when threshold is exceeded
%     TIME      
%         .it             : timestep index (starts at it=0 at model start)
%     i_mc                : index of the active coastline element       
% 
% OUTPUT:
%     COAST
%         .x              : x-coordinates of considered coastal section [m]
%         .y              : y-coordinates of considered coastal section [m]
%         .n              : number of points of the considered coastal section
%         .ds             : grid cell size [m]
%         .xq             : x-coordinates of transport points in section (m)
%         .yq             : y-coordinates of transport points in section (m)
%         .nq             : number of transport points in the considered coastal section.
%         .dsq            : grid cell size for QS-points [m]
%         .x_mc           : x-coordinate of coastline (all sections)
%         .y_mc           : y-coordinate of coastline (all sections)
%         .n_mc           : number of coastal sections
%         .cyclic         : index describing whether the considered coastline section is cyclic
%         .clockwise      : index describing whether the considered coastline section is cyclic
%         .idgrid         : identifiers storing the index of the gridcells in the former administration for the new grid (e.g. showing a number with a 0.5 is a cell was doubled)
%         .idgridnr       : identifiers storing cells that were split or doubled (e.g. showing 0 if a cell is removed and 2 if it is doubled)
%         .gridchange     : status indicator showing whether the grid has changed in this timestep-iteration
% 
% Example :
%         COAST.x_mc=[0,40,1000,1200,2500,4022,4023];
%         COAST.y_mc=[100,101,100,100,200,300,320];
%         TIME.it=10;
%         i_mc=1;
%         [COAST]=make_sgrid_mc(COAST,TIME,i_mc);
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

    eps=0.1;
    sqrt2=sqrt(2.);
    %% find number of sections
    [x,y]=get_one_polygon( COAST.x_mc,COAST.y_mc,i_mc );  
    clockwise=get_clockpoly(x,y);
    COAST.gridchange=0;

    %% METHOD 1: preserve existing points as much as possible by splitting cells -> a quick method wherein only relevant cells are changed with low diffusion, but with a less smooth grid
    if COAST.griddingmethod==1
        
        % iteratively make the grid
        niter=1;
        if TIME.it==0 
            niter=10; % at t=0 it needs more iterations
            % make sure to update the grid interpolation at t0
            COAST.gridchange=1;
        end 
        
        for mm=1:niter
            idgridding=[];
            % interpolate grid cell size
            [ds0i,ds0imean]=get_gridsize(x,y,COAST.ds0);
            % create grid
            s=[0,cumsum(hypot(diff(x),diff(y)))];
            ds=diff(s);
            if ~isempty(sum(ds==0))  
                idu=setdiff([1:length(s)],find(ds==0));
                x=x(idu);
                y=y(idu);
                s=s(idu);
                ds=ds(idu(1:end-1));
                [ds0i,ds0imean]=get_gridsize(x,y,COAST.ds0);
            end
            
            if max(ds./ds0imean>2)==1 || max(ds./ds0imean<0.5)==1
            %if max(ds)>2*COAST.ds0 || min(ds)<COAST.ds0/2
                COAST.gridchange=1;
                % compute number of grid cells that can fit in-between each xy-point.
                sn0=ds./ds0imean;

                % make sure last x,y point is accounted for by adding the relative fraction of a grid cell of adjacent cells (near the end, until it is at least > 0.5*ds0)
                if sn0(end)<0.5
                    sumsn0=cumsum(fliplr(sn0));
                    sn0half=find(sumsn0>0.5,1);
                    if ~isempty(sn0half)
                        % make sure to end the coast at the end point and not 1 or more grid points earlier
                        % to do so the gridcell weight (or fraction) of the cells before it are attributed to the last point. 
                        sn0(length(sn0)-sn0half+1:end)=0;
                        sn0(end)=sumsn0(sn0half);
                    else
                        % always make at least 2 points! 
                        % even if it was just one
                        sn0(end)=0.5;
                    end
                end
                
                % make the new distance table by each time adding the number of grid cells to the new xy-point
                % if this number of grid cells (sn0) is 0 then no grid cells are added, and the xy-point is omitted.
                sn0rest=0;
                snew=s(1);            
                for ii=1:length(sn0)
                    s2=[];
                    id2=[];
                    idgrid=1;
                    sn=round(sn0(ii)+sn0rest);
                    idgridnr(1,ii)=sn;
                    sn0rest=max(sn0(ii)+sn0rest-sn,0);
                    if sn>0
                        s2=[snew(end):(s(ii+1)-snew(end))/sn:s(ii+1)];   
                        id2=idgrid(end)+[1:1:sn]/sn;
                    end
                    snew=[snew,s2(2:end)];
                    idgrid=[idgrid,id2];
                end
                x=interp1(s,x,snew);
                y=interp1(s,y,snew);
                s=snew;
            else
                idgrid=[1:length(x)];
                idgridnr=ones(1,length(x));
            end
            %figure;plot(x,y,'k-','LineWidth',2);hold on;plot(x,y,'ks','LineWidth',2);plot(x,ynew,'r.-','LineWidth',2);hold on;
        end
        
        % identifiers storing the index of the gridcells in the former administration for the new grid (e.g. showing a number with a 0.5 is a cell was doubled)
        COAST.idgrid=idgrid;
        % identifiers storing cells that were split or doubled (e.g. showing 0 if a cell is removed and 2 if it is doubled)
        COAST.idgridnr=idgridnr;
        
    %% METHOD 2: regenerate the grid every time the threshold is exceeded -> results in a smoother grid, but is slower and generates more diffusion
    elseif COAST.griddingmethod==2
        ds=hypot(diff(x),diff(y));
        s=[0,cumsum(ds)];
        ds0=COAST.ds0;
        if ~isscalar(ds0)
            ds0=min(ds0(:,3));
        end
        
        %if max(ds./ds0imean>2)==1 || TIME.it==0 
        if max(ds)>2*ds0 || TIME.it==0 
            COAST.gridchange=1;
            % remove double points
            IDunique = [true,ds>eps];
            x0=x(IDunique);
            y0=y(IDunique);
            s0=s(IDunique);
            % compute distance along line
            ns=ceil(s0(end)/ds0);
            ds1=s0(end)/ns;
            s=[0:ds1:s0(end)];
            % interpolate x-values
            x=interp1(s0,x0,s);
            y=interp1(s0,y0,s);
        end
        snew=s;
        i=2;
        %[ds0i,ds0imean]=get_gridsize(x,y,COAST.ds0);
        while i<=length(snew)            
            ds2=snew(i)-snew(i-1);
            if ds2<ds0/sqrt2 && snew(end)>=ds0/sqrt2
                COAST.gridchange=1;
                %throw out point i
                if i>2 && i<length(snew)
                    snew=[snew(1:i-1),snew(i+1:end)];
                elseif i==2 && ds2<ds0/4 && length(snew)>2
                    snew=[snew(1),snew(3:end)];
                elseif ds2<ds0/4
                    snew=[snew(1:i-2),snew(end)];
                else
                    i=i+1;
                end
            elseif ds2>ds0*sqrt2 
                COAST.gridchange=1;
                %insert point i                
                snew=[snew(1:i-1),.5*(snew(i-1)+snew(i)),snew(i:end)];
                i=i+1;
            else
                i=i+1;
            end
        end
        
        snew(2:end-1)=COAST.smoothfac*snew(1:end-2)+(1.-2*COAST.smoothfac)*snew(2:end-1)+COAST.smoothfac*snew(3:end);
        if length(snew)~=length(s) || sum(s)~=sum(snew)
            if clockwise==1 || min(s)<ds0/5
            x=interp1(s,x,snew);
            y=interp1(s,y,snew);
            s=snew;
            end
        end
    end
    
    %% Check if coastal section is cyclic
    cyclic=get_cyclic(x,y,COAST.ds0);
    if cyclic
        x(end)=x(1);
        y(end)=y(1);
    end
    
    %% Compute grid cell size
    ds0=diff(s);
    if ~isempty(ds0)
        if cyclic
            ds=[(ds0(1)+ds0(end))/2,(ds0(1:end-1)+ds0(2:end))/2,(ds0(1)+ds0(end))/2]; 
        else
            ds=[ds0(1),(ds0(1:end-1)+ds0(2:end))/2,ds0(end)]; 
        end
        
        %% Check if coastal section is clockwise (clockwise=1 means a coastline, while clockwise=0 means a lake)
        clockwise=get_clockpoly(x,y);
        if ~cyclic
            clockwise=1;
        end
        
        % compute locations of transport points
        % and add extra boundary point
        xq=(x(1:end-1)+x(2:end))/2;
        yq=(y(1:end-1)+y(2:end))/2;
        if cyclic
            xq=[xq(end),xq,xq(1)];
            yq=[yq(end),yq,yq(1)];
        else
            if COAST.BNDgroyne(i_mc,1)==0       % extend grid when it is an open coast without a structure
                xq=[1.5*x(1)-0.5*x(2),xq];
                yq=[1.5*y(1)-0.5*y(2),yq];
            else                                % use first coastline point in case a groyne is present
                xq=[x(1),xq];
                yq=[y(1),yq];
            end
            if COAST.BNDgroyne(i_mc,2)==0       % extend grid when it is an open coast without a structure
                xq=[xq,1.5*x(end)-0.5*x(end-1)];
                yq=[yq,1.5*y(end)-0.5*y(end-1)];
            else                                % use last coastline point in case a groyne is present
                xq=[xq,x(end)];
                yq=[yq,y(end)];
            end
        end
    else
        ds=[];
        xq=[];
        yq=[];
        clockwise=0;
        cyclic=0;
    end
    
    % update structures
    COAST.i_mc=i_mc;
    COAST.x=x;
    COAST.y=y;
    COAST.n=length(x);
    COAST.s=s;
    COAST.ds=ds;
    COAST.xq=xq;
    COAST.yq=yq;
    COAST.dsq=[ds(1),(ds(1:end-1)+ds(2:end))/2,ds(end)];
    COAST.nq=length(xq);
    COAST.cyclic=cyclic;
    COAST.clockwise=clockwise;
    
    %% insert x,y back into COAST.x_mc,COAST.y_mc
    [COAST.x_mc,COAST.y_mc]=insert_section(COAST.x,COAST.y,COAST.x_mc,COAST.y_mc,COAST.i_mc);
    COAST.n_mc=length(find(isnan(COAST.x_mc)))+1;
end


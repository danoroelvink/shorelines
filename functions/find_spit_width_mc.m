function [ COAST,SPIT ] = find_spit_width_mc(i_mc,COAST,WAVE,SPIT,STRUC,TRANSP,TIME)
% function [ COAST,SPIT ] = find_spit_width_mc(i_mc,COAST,WAVE,SPIT,STRUC,TRANSP)
%
% <Summary of this function goes here>
% For each section
%    Select coastline points not in shadow of coastline or structures
%    For each selected coastline point
%        Construct landward normal
%        For each coastline section
%            Find intersections with landward normal
%            Compute distance and save distance, point numbers, weights and section number
%            Select shortest distance 
%            Compute overwash displacement at coastline point; 
%               Simple function based on width or function of wave condition, water level
%               check cyclic start or end point
%            Compute displacement backside points; 
%               Spread over adjacent points 
%               check cyclic start or end point

%
% INPUT:
%      COAST.x_mc       x-coordinates of grid
%
% OUTPUT:
%      COAST.x_mc       x-coordinates of grid (corrected for overwash)
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
    
    %% INITIALIZE VARIABLES
    nrpoints_inbetween = 2;     % minimum number of points in-between overwash points 
    eps=5;
    yesplot=0;
    SPIT.width = zeros(size(COAST.x_mc));
    SPIT.spit  = false(size(COAST.x_mc));
    shadow     = false(size(COAST.x_mc));
    shadow_h   = false(size(COAST.x_mc));
    dxt        = zeros(size(COAST.x_mc));
    dyt        = zeros(size(COAST.x_mc));
    lenw       = 10000;
    n          = length(COAST.x_mc);

    %% find out number of sections n_mc and start and end indices of all sections i1 and i2
    i1=zeros(COAST.n_mc,1);
    i2=zeros(COAST.n_mc,1);
    
    for ii=1:COAST.n_mc
        [ ~,~,~,i1(ii),i2(ii) ] = get_one_polygon( COAST.x_mc,COAST.y_mc,ii);
    end
    
    %% In the following, ii denotes the section from where we check the spit width
    %  and jj the section of the point on the other side of the spit

    [COAST,WAVE,~]=get_segmentdata(COAST,WAVE,TRANSP,i_mc);
    if length(COAST.x)<3
        return
    end
    % COAST.cyclic has to be recomputed
    xc=COAST.x_mc(i1(i_mc):i2(i_mc));
    yc=COAST.y_mc(i1(i_mc):i2(i_mc));
    COAST.cyclic=get_cyclic(xc,yc,COAST.ds0);
    if 0 %COAST.cyclic
        COAST.x_mc(i2(i_mc))=COAST.x_mc(i1(i_mc));
        COAST.y_mc(i2(i_mc))=COAST.y_mc(i1(i_mc));
    end 
    
    % compute PHIWI2 at nodes from PHIWI at vertices.
    phiwi2=mod(atan2d(sind(0.5*(WAVE.PHI(1:end-1)+WAVE.PHI(2:end))),cosd(0.5*(WAVE.PHI(1:end-1)+WAVE.PHI(2:end)))),360.0);
    if ~isempty(phiwi2)
        phiwi2=[phiwi2(1),phiwi2,phiwi2(end)];
        if COAST.cyclic
            phiwi2(1)=mod(atan2d(sind(0.5*(phiwi2(1)+phiwi2(end))),cosd(0.5*(phiwi2(1)+phiwi2(end)))),360.0);
            phiwi2(end)=phiwi2(1);
        end 
    end
    
    for i=i1(i_mc):i2(i_mc)
        iloc=i-i1(i_mc)+1;
        
%         if ~isnan(COAST.x_mc(i))  % if the point not Nan -> this can
%         never happen?
            if COAST.cyclic 
                im1=i-1;if im1<i1(i_mc);im1=i2(i_mc)-1;end
                ip1=i+1;if ip1>i2(i_mc);ip1=i1(i_mc)+1;end
            else
                im1=max(i-1,i1(i_mc));
                ip1=min(i+1,i2(i_mc));
            end
            
            % construct line in wave direction and check for shadowing
            phic=360-atan2d(COAST.y_mc(ip1)-COAST.y_mc(im1),COAST.x_mc(ip1)-COAST.x_mc(im1));
            xl=[COAST.x_mc(i)+eps*sind(phic),COAST.x_mc(i)-SPIT.spit_width*sind(phic)];  
            yl=[COAST.y_mc(i)+eps*cosd(phic),COAST.y_mc(i)-SPIT.spit_width*cosd(phic)];  
            
            % check for shadowing by other section
            try
               xw=[COAST.x_mc(i)+eps*sind(phiwi2(iloc)),COAST.x_mc(i)+lenw*sind(phiwi2(iloc))]; 
               yw=[COAST.y_mc(i)+eps*cosd(phiwi2(iloc)),COAST.y_mc(i)+lenw*cosd(phiwi2(iloc))]; 
            end
            
            % Check for shadowing
            [xx1,yy1]=get_intersections(COAST.x_mc,COAST.y_mc,xl,yl);
            if isempty(xx1), continue, end
            if xx1==-1e10, continue, end
            P1=[xx1;yy1];
            
            if (isoctave)
               P1=get_uniquetoloct(P1',1d-12,'rows');
               P1=P1(:)'; 
            else    
               P1=uniquetol(P1','byrows',true)';
               P1=P1';
            end
            [xx2,~]=get_intersections(COAST.x_mc,COAST.y_mc,xw,yw);
            shadow(i)=~isempty(xx2); 
            if ~isempty(STRUC.x_hard)
                [xx3,~]=get_intersections(STRUC.x_hard,STRUC.y_hard,xw,yw);
                [xx4,~]=get_intersections(STRUC.x_hard,STRUC.y_hard,xl,yl);
                shadow_h(i)=~isempty(xx3);  % shadow by hard structure?
                shadow_sp(i)=~isempty(xx4); % structure in spit?
            else
                shadow_h(i)=0;
                shadow_sp(i)=0;
            end
            
            % Set spit overwash point to 0 if there is some shadowing present in considered cell
            SPIT.spit(i)=length(xx1)==2 && ~shadow(i) && ~shadow_h(i) && ~shadow_sp(i);  % if the spit width > critcal & not shadowed
            if SPIT.spit(i)
                SPIT.width(i)=hypot(xx1(1)-xx1(2),yy1(1)-yy1(2));   % if so calculate the actual width
            end
            if yesplot
                figure(10);
                if size(P1,2)*size(P1,1)>2 && size(P1,2)>1
                    plot(COAST.x_mc,COAST.y_mc,COAST.x_mc(i),COAST.y_mc(i),'^g',xl,yl,'ro-',xx1(1),yy1(1),'*k',xx1(2),yy1(2),'*b',xw,yw,'mo-');
                    title(['width=',num2str(SPIT.width(i))]);
                end
            end
                    
        %% 
        if SPIT.width(i)>0.1 && SPIT.width(i)<SPIT.spit_width ...
            && length(xx1)>1 && shadow(i)==0 && shadow(i)==0  % if applicable %here we miss the shadowing by structures
            
            %% to find the second closest point on the other side
            if hypot(COAST.x_mc(i)-xx1(1),COAST.y_mc(i)-yy1(1))>hypot(COAST.x_mc(i)-xx1(2),COAST.y_mc(i)-yy1(2))  %determine the closest point on the other side
                P_diff=(hypot(COAST.x_mc-xx1(1),COAST.y_mc-yy1(1)));
            else
                P_diff=(hypot(COAST.x_mc-xx1(2),COAST.y_mc-yy1(2)));
            end
            [P_dist,P_idx]=sort(P_diff); %sort them
            ip=P_idx(1); % the closest
            ipp=P_idx(2); % the second closest 
            mindist=P_dist(1);
            mindist2=P_dist(2);
            al=1-(mindist/(mindist+mindist2));
            
            % determine the point before and the point after
            if ip>0&&ip<=n && abs(ip-i)>nrpoints_inbetween
                if COAST.cyclic
                    im1=i-1;if im1<i1(i_mc);im1=i2(i_mc)-1;end
                    ip1=i+1;if ip1>i2(i_mc);ip1=i1(i_mc)+1;end
                else
                    im1=max(i-1,i1(i_mc));
                    ip1=min(i+1,i2(i_mc));
                end
                if yesplot
                    figure(10);
                    plot(COAST.x_mc,COAST.y_mc,xl,yl,'.-',xx1,yy1,'ro',xx1(2),yy1(2),'k+','linewidth',2);
                    axis equal
                    title(['width ',num2str(round(SPIT.width(i)))]);
                    hold on
                    plot(COAST.x_mc(ip),COAST.y_mc(ip),'ok',COAST.x_mc(im1),COAST.y_mc(im1),'ob',COAST.x_mc(ip1),COAST.y_mc(ip1),'ob')
                end
                
                %% shift seaward side of spit
                if ~(COAST.cyclic&&i==i2(i_mc))
                    phic=360-atan2d(COAST.y_mc(ip1)-COAST.y_mc(im1),COAST.x_mc(ip1)-COAST.x_mc(im1));
                    try
                       if SPIT.OWtimescale>0
                          dn=-TIME.dt/SPIT.OWtimescale*(SPIT.spit_width-SPIT.width(i))*cosd(phiwi2(iloc)-phic)*(SPIT.Dbb+SPIT.Bheight)/(SPIT.Dsf+SPIT.Bheight); % <- sometimes issue arise here : "Index exceeds the number of array elements (3)."
                       else
                          dn=-SPIT.OWscale*(SPIT.spit_width-SPIT.width(i))*cosd(phiwi2(iloc)-phic)*(SPIT.Dbb+SPIT.Bheight)/(SPIT.Dsf+SPIT.Bheight); % <- sometimes issue arise here : "Index exceeds the number of array elements (3)."
                       end
                    catch
                       disp('gotcha fella') 
                    end
                    ds=hypot(COAST.x_mc(ip1)-COAST.x_mc(im1),COAST.y_mc(ip1)-COAST.y_mc(im1));
                    % make sure element gridsize is not too small! Otherwise do not use the spit width function. Element will be adjusted later in cleanup_nans routine anyways.
                    if ds<eps
                        dn=0;ds=1;
                    end
                    dx=-dn*(COAST.y_mc(ip1)-COAST.y_mc(im1))/hypot(COAST.x_mc(ip1)-COAST.x_mc(im1),COAST.y_mc(ip1)-COAST.y_mc(im1));
                    dy= dn*(COAST.x_mc(ip1)-COAST.x_mc(im1))/hypot(COAST.x_mc(ip1)-COAST.x_mc(im1),COAST.y_mc(ip1)-COAST.y_mc(im1));
                    dA=dn*ds;
                    dxt(i)=dxt(i)+dx;
                    dyt(i)=dyt(i)+dy;
                    if COAST.cyclic&&i==i1(i_mc) %to keep the cyclic section closed
                        COAST.x_mc(i2(i_mc))=COAST.x_mc(i);
                        COAST.y_mc(i2(i_mc))=COAST.y_mc(i);
                        dxt(i2(i_mc))=dxt(i);
                        dyt(i2(i_mc))=dyt(i);
                    end
                    %% shift landward side of spit
                    % what is jj?
                    for j=1:COAST.n_mc   % determine the section where the landward point exist
                        if ip>=i1(j) && ip<=i2(j)
                            jj=j;
                            break
                        end
                    end
                    % is jj cyclic?
                    cyclic_jj= hypot(COAST.x_mc(i2(jj))-COAST.x_mc(i1(jj)),COAST.y_mc(i1(jj))-COAST.y_mc(i2(jj)))<ds;
                    if cyclic_jj;COAST.x_mc(i2(jj))=COAST.x_mc(i1(jj));COAST.y_mc(i2(jj))=COAST.y_mc(i1(jj));end
                    if cyclic_jj
                        im1=ip-1;if im1<i1(jj);im1=i2(jj)-1;end
                        ip1=ip+1;if ip1>i2(jj);ip1=i1(jj)+1;end
                    else
                        im1=max(ip-1,i1(jj));
                        ip1=min(ip+1,i2(jj));
                    end
                    ds=hypot(COAST.x_mc(ip1)-COAST.x_mc(im1),COAST.y_mc(ip1)-COAST.y_mc(im1));
                    dn=dA/ds*((SPIT.Dsf+SPIT.Bheight)/(SPIT.Dbb+SPIT.Bheight))^2;
                    % make sure element gridsize is not too small! Otherwise do not use the spit width function. Element will be adjusted later in cleanup_nans routine anyways.
                    if ds<eps
                        dn=0;ds=1;
                    end
                    dx= dn*(COAST.y_mc(ip1)-COAST.y_mc(im1))/ds;
                    dy=-dn*(COAST.x_mc(ip1)-COAST.x_mc(im1))/ds;
                    dxt(ip)=dxt(ip)+dx*al;
                    dyt(ip)=dyt(ip)+dy*al;
                    dxt(ipp)=dxt(ipp)+dx*(1-al);
                    dyt(ipp)=dyt(ipp)+dy*(1-al);
                    
                    if cyclic_jj
                        if ip==i1(jj)
                            COAST.x_mc(i2(jj))=COAST.x_mc(ip);
                            COAST.y_mc(i2(jj))=COAST.y_mc(ip);
                            dxt(i2(jj))=dxt(ip);
                            dyt(i2(jj))=dyt(ip);
                        elseif ip==i2(jj)
                            COAST.x_mc(i1(jj))=COAST.x_mc(ip);
                            COAST.y_mc(i1(jj))=COAST.y_mc(ip);
                            dxt(i1(jj))=dxt(ip);
                            dyt(i1(jj))=dyt(ip);
                        end
                    end
                end
                if yesplot
                    plot(COAST.x_mc,COAST.y_mc,'--k')
                    hold off
                    drawnow
                    %%pause
                end
            end
        end
    end
    
    COAST.x_mc=COAST.x_mc+dxt;
    COAST.y_mc=COAST.y_mc+dyt;
    %% extract x and y
    COAST.x = COAST.x_mc(i1(i_mc):i2(i_mc));
    COAST.y = COAST.y_mc(i1(i_mc):i2(i_mc));
    
end

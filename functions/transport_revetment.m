function [TRANSP,COAST]=transport_revetment(COAST,STRUC,TRANSP,TIME,WAVE)
%function [TRANSP,COAST]=transport_revetment(COAST,STRUC,TRANSP,TIME,WAVE)
%
% Reduces transports if the coastline is near a revetment; limits
% the transport out of a cell to the volume of the cell divided by the 
% timestep.
%
% INPUT: 
%   COAST
%         .x           : x-coordinate of coastline [m] (only current section)
%         .y           : y-coordinate of coastline [m] (only current section)
%         .s           : alongshore distance along the coastline [m]
%         .ds0         : default grid cell size [m]
%         .ds          : grid cell size [m]
%         .h0          : active height of the profile [m]
%   STRUC
%         .xrevet      : x-coordinate of revetments [m]
%         .yrevet      : y-coordinate of revetments [m]
%         .iterrev     : iterations used for determining transports along revetments [-]
%   TRANSP
%         .QS          : transport rates [m3/yr]
%         .critwidth   : critical beach width below which transport is reduced [m]
%         .xsedlim     : x-coordinate of regions with sediment limitation [m]
%         .ysedlim     : y-coordinate of regions with sediment limitation [m]
%         .widthsedlim : critical width at regions with sediment limitation, below which transport is reduced [m]
%   TIME
%         .it          : timestep index since model start [-]
%         .adt         : time step [year]    
%
% OUTPUT:
%   TRANSP
%         .QS          : modified transport
%         .idrev       : index with locations of revetments
%   COAST
%         .x           : updated x-coordinate of coastline (only current section)
%         .y           : updated y-coordinate of coastline (only current section)
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2022 IHE Delft & Deltares
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

    %% Initializations
    eps=1e-6;
    nw       = size(TRANSP.QS,1); % number of wave conditions taken along at this timestep (can be more than 1 in case of simultaneous wave conditions)
    nq       = size(TRANSP.QS,2); % number of coastline points
    nc       = length(COAST.x);    % number of coastline points
    x        = COAST.x;                
    y        = COAST.y;               
    s        = COAST.s; 
    ds0      = COAST.ds0;
    ds       = COAST.ds; % cell size
    QS       = TRANSP.QS;          
    iterrev  = STRUC.iterrev;
    h0       = COAST.h0;
    reducfac = ones(1,nq);
    dnrev    = ones(1,nc)*1e4;     % waterline distance in front of revetment (very wide by default)
    dnrevq   = ones(1,nq)*1e4;     % waterline distance in front of revetment on transport grid cells (very wide by default)
    TRANSP.idrev=zeros(1,nq);      % QS-indices with a revetment -> used to communicate to 'upwind correction' and 'bypass function'
    if length(h0)==1
        h0=repmat(h0,[1,nc]);
    end
    if ~isscalar(ds0)
        ds0=min(ds0(:,3));
    end

    %% get number of revetments
    [x_str,y_str,n_str]=get_one_polygon(STRUC.xrevet,STRUC.yrevet,1);
    [x_sed,y_sed,n_sed]=get_one_polygon(TRANSP.xsedlim,TRANSP.ysedlim,1);
    
    %% for all revetments
    TRANSP.reducfac=ones(1,nq);
    TRANSP.idrevq=zeros(1,nq);
    if n_str>0 || n_sed>0
        for i_str = 1:n_str+n_sed
            if i_str<=n_str
                % get the revetment coordinates
                [x_rev,y_rev,n_str]=get_one_polygon(STRUC.xrevet,STRUC.yrevet,i_str);
                wcrit = TRANSP.critwidth;
                wreduce = 0;
            else
                % get 'virtual revetment' with sediment limitation
                i_sed=i_str-n_str;
                [x_sed,y_sed,n_sed]=get_one_polygon(TRANSP.xsedlim,TRANSP.ysedlim,i_sed);
                wcrit = mean(get_one_polygon(TRANSP.widthsedlim,i_sed));
                wreduce = wcrit;
                x_rev=x_sed;
                y_rev=y_sed;
            end
            
            % compute distances from the coast to the revetment
            usepolydistance=1;
            if usepolydistance==0
                [dmin,x_dmin,y_dmin,is_vertex,idc]=p_poly_dist(x,y,x_rev,y_rev);
                % Convert from column to row vectors
                dmin=dmin';x_dmin=x_dmin';y_dmin=y_dmin';is_vertex=is_vertex';idc=idc';
                % determine coast points within the revetment area
                inside=~(idc==0&is_vertex|idc==length(x_rev)-1&is_vertex);
                a=[dif2(x);dif2(y)];             % vector in direction of coastline
                b=[x_dmin-x;y_dmin-y];           % vector from coastline point to closest
                                                 % point on revetment
                % cross product: if negative then coast seaward of revetment
                cr=a(1,:).*b(2,:)-a(2,:).*b(1,:);
                % make dmin positive when coast seward of revetment
                dmin(inside)=-dmin(inside).*sign(cr(inside));
            else
                % extend with 0% of a standard grid cell size to get some extra coverage (avoid undermining)
                dx=diff(x_rev);dx=[dx(1),dx(end)];
                dy=diff(y_rev);dy=[dy(1),dy(end)];
                L=0.05*ds0./((dx.^2+dy.^2).^0.5);  % a value of 0 is used, which means the revetment is not extended
                x_revext=[x_rev(1)-dx(1)*L(1),x_rev,x_rev(end)+dx(end)*L(end)];
                y_revext=[y_rev(1)-dy(1)*L(1),y_rev,y_rev(end)+dy(end)*L(end)];
                
                % compute distance to (extended) revetment
                [dmin,x_dmin,y_dmin]=get_polydistance(x,y,x_revext,y_revext);
                inside=~isnan(dmin);
                dmin=-dmin;
            end
            
            % computing the cross-shore distance to the revetment for the coastline points (dnrev) and for the qs-points (dnrevq)
            dnrev(inside)=max(dmin(inside),0);
            dnrevq=[dnrev(1),(dnrev(1:end-1)+dnrev(2:end))/2,dnrev(end)];
            
            % coastline behind revetment
            insideq=[inside;0] | [0;inside];
            behind=inside&dmin<=wreduce; 
            behindq=[behind;0] | [0;behind]; 
            
            % at t0 the coast is fixed to revetment
            if sum(behind)>0 && TIME.it==0
                x(behind)=x_dmin(behind);
                y(behind)=y_dmin(behind);
                dnrevq(insideq)=max(dnrevq(insideq),0);
                %% insert x,y back into COAST.x_mc,COAST.y_mc
                COAST.x=x;
                COAST.y=y;
                [COAST.x_mc,COAST.y_mc]=insert_section(COAST.x,COAST.y,COAST.x_mc,COAST.y_mc,COAST.i_mc);
            end
            % reduction factor
            reducfac(insideq)=max(min(dnrevq(insideq)/wcrit,1),0);               
            TRANSP.idrev=TRANSP.idrev | behindq(:)';    % QS-indices with a revetment -> used to communicate to upwind correction and bypass function
            TRANSP.reducfac=reducfac;
            % debug plot
            % figure;plot(COAST.x,COAST.y,'k.-');hold on;plot(COAST.xq(insideq),COAST.yq(insideq),'gs');hold on;plot(COAST.x(inside),COAST.y(inside),'r+');
            
            % reduce transports to avoid emptying coastal cells in front of revetments
            iddivergence=zeros(1,nq-1);
            for kk=1:nw
            iddivergence=iddivergence | (QS(kk,2:end)>1e-6 & QS(kk,1:end-1)<1e-6);
            end
            dnrev(iddivergence)=dnrev(iddivergence)/2;
            
            % reduce transports to mimic limitation to bypass
            for kk=1:nw
                QS(kk,:)=QS(kk,:).*reducfac;
                if isempty(WAVE.Prob)
                    Prob=1;
                else
                    Prob=WAVE.Prob(kk)/sum(WAVE.Prob(1:nw));
                end
                for iter=1:iterrev 
                    QS0=QS(kk,:);
                    for i=2:nq
                        if QS(kk,i)>=0
                            QS(kk,i)=min(QS(kk,i),QS0(i-1)+dnrev(i-1)*ds(i-1)*h0(i-1)/TIME.adt/Prob);
                            %QS(kk,i)=min(QS(kk,i),dnrev(i-1)*ds(i-1)*h0(i-1)/TIME.adt);        % <- basic approach, also works well, but bypass transport is a linear function of the time step
                        end
                    end
                    for i=[nq-1:-1:1]
                        if QS(kk,i)<0
                            QS(kk,i)=max(QS(kk,i),QS0(i+1)-dnrev(i)*ds(i)*h0(i)/TIME.adt/Prob);
                            %QS(kk,i)=max(QS(kk,i),-dnrev(i)*ds(i)*h0(i)/TIME.adt);        % <- basic approach, also works well, but bypass transport is a linear function of the time step
                        end
                    end
                end
                
                %% shadowing of the revetment
                % TRANSP.QS(dnrevq==1e4 & TRANSP.shadowS_rev)=0;
                TRANSP.idrevq=(dnrevq<100);
                TRANSP.idrevs=(dnrev<100);
                % idrevq2=[0,cumsum(abs(diff(dnrevq<1e4)))];
                idshadowrevneg=intersect(find((TRANSP.shadowS_revq(kk,:) | [TRANSP.shadowS_rev(kk,:),0]) ),[1:find(TRANSP.idrevq==1,1)]);
                idshadowrevpos=intersect(find((TRANSP.shadowS_revq(kk,:) | [0,TRANSP.shadowS_rev(kk,:)]) ),[find(TRANSP.idrevq==1,1,'last'):nq]);
                % idshadowrevpos=find(TRANSP.shadowS_rev & TRANSP.QS>=0 & idrevq2==0);
                if length(idshadowrevneg)>=2
                    QS(kk,idshadowrevneg)=min(QS(kk,min(idshadowrevneg(end)+1,nq)),0)*[0:1/(length(idshadowrevneg)-1):1];
                end
                if length(idshadowrevpos)>=2
                    QS(kk,idshadowrevpos)=max(QS(kk,max(idshadowrevpos(1)-1,1)),0).*[1:-1/(length(idshadowrevpos)-1):0];
                end
            end
        end
        
        COAST.x   = x;
        COAST.y   = y;
        TRANSP.QS  = QS;
    end
end

function dif=dif2(x)
    dif=zeros(size(x));
    dif(1)=x(2)-x(1);
    dif(2:end-1)=(x(3:end)-x(1:end-2))/2;
    dif(end)=x(end)-x(end-1);
end

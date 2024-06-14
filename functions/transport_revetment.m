function [TRANSP,COAST]=transport_revetment(COAST,STRUC,TRANSP,TIME)
%function [TRANSP,COAST]=transport_revetment(COAST,STRUC,TRANSP,TIME)
%
% Reduces transports if the coastline is near a revetment; limits
% the transport out of a cell to the volume of the cell divided by the 
% timestep.
%
% INPUT :
%   COAST
%         .x         : x-coordinate of coastline (only current section)
%         .y         : y-coordinate of coastline (only current section)
%   STRUC
%         .x_revet   : x-coordinate of revetments
%         .y_revet   : y-coordinate of revetments
%   TRANSP
%         .xS        : x-coordinate of QS-points
%         .yS        : y-coordinate of QS-points
%         .QS        : transport
%         .crit_width: critical beach width below which transport is reduced
%   TIME
%         .dt        : latest time step
%
% OUTPUT:
%   TRANSP
%         .QS        : modified transport
%   COAST
%         .x         : updated x-coordinate of coastline (only current section)
%         .y         : updated y-coordinate of coastline (only current section)
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
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses
%   --------------------------------------------------------------------
%% Initializations
    eps=1e-6;
    x        = COAST.x;                
    y        = COAST.y;               
    s        = COAST.s; 
    ds0      = COAST.ds0;
    ds       = COAST.ds; % cell size
    QS       = TRANSP.QS;              
    wcrit    = TRANSP.crit_width;   
    iterrev  = STRUC.iterrev;
    h0       = COAST.h0;
    reducfac = ones(size(QS));
    dnrev    = ones(size(x))*1e4;       % waterline distance in front of revetment (very wide by default)
    dnrevq   = ones(size(QS))*1e4;   % waterline distance in front of revetment (very wide by default)
    TRANSP.idrev=zeros(1,length(QS));    % QS-indices with a revetment -> used to communicate to 'upwind correction' and 'bypass function'
    if length(h0)==1
        h0=repmat(h0,size(x));
    end
    if ~isscalar(ds0)
        ds0=min(ds0(:,3));
    end

    % get number of revetments
    [x_str,y_str,n_str]=get_one_polygon(STRUC.x_revet,STRUC.y_revet,1);
    [x_sed,y_sed,n_sed]=get_one_polygon(TRANSP.x_sedlim,TRANSP.y_sedlim,1);
    
    %% for all revetments
    if n_str>0 || n_sed>0
        for i_str = 1:n_str+n_sed
            if i_str<=n_str
                % get the revetment
                [x_rev,y_rev,n_str]=get_one_polygon(STRUC.x_revet,STRUC.y_revet,i_str);
                wcrit = TRANSP.crit_width;
                wreduce = 0;
            else
                % get 'virtual revetment' with sediment limitation
                i_sed=i_str-n_str;
                [x_sed,y_sed,n_sed]=get_one_polygon(TRANSP.x_sedlim,TRANSP.y_sedlim,i_sed);
                wcrit = mean(get_one_polygon(TRANSP.width_sedlim,i_sed));
                wreduce = wcrit;
                x_rev=x_sed;
                y_rev=y_sed;
            end
            
            % compute distances from the coast to the revetment
            usepolydistance=1;
            if usepolydistance==0
                [dmin,x_dmin, y_dmin, is_vertex, idc]=p_poly_dist(x,y,x_rev,y_rev);
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
                dnrev(inside)=max(dmin(inside),0);

                % coastline behind revetment; is fixed to revetment
                behind=inside&dmin<wreduce;            
                if sum(behind)>0
                x(behind)=x_dmin(behind);
                y(behind)=y_dmin(behind);
                end
                %% Reduction factor
                reducfac(inside)=max(min(dnrev(inside)/wcrit,1),0);
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
                dnrev(inside)=max(dmin(inside),0);
                dnrevq=[dnrev(1),(dnrev(1:end-1)+dnrev(2:end))/2,dnrev(end)];
                insideq=[inside;0] | [0;inside];
                
                % coastline behind revetment; is fixed to revetment
                behind=inside&dmin<=wreduce; 
                behindq=[behind;0] | [0;behind]; 
                if sum(behind)>0 && TIME.it==0
                    x(behind)=x_dmin(behind);
                    y(behind)=y_dmin(behind);
                    dnrevq(insideq)=max(dnrevq(insideq),0);
                end
                reducfac(insideq)=max(min(dnrevq(insideq)/wcrit,1),0);
                
                TRANSP.idrev=TRANSP.idrev | behindq(:)';    % QS-indices with a revetment -> used to communicate to 'upwind correction' and 'bypass function
            end
            
            %% Debug plot
            % figure;plot(COAST.x,COAST.y,'k.-');hold on;plot(COAST.xq(insideq),COAST.yq(insideq),'gs');hold on;plot(COAST.x(inside),COAST.y(inside),'r+');
            
            %% Reduce transports to mimic limitation to bypass
            QS=QS.*reducfac;
            
            %% Reduce transports to avoid emptying coastal cells in front of revetments
            iddivergence=QS(2:end)>0 & QS(1:end-1)<0;
            dnrev(iddivergence)=dnrev(iddivergence)/2;
            for iter=1:iterrev 
                QS0=QS;
                for i=2:length(QS)
                    if QS(i)>=0
                        QS(i)=min(QS(i),QS0(i-1)+dnrev(i-1)*ds(i-1)*h0(i-1)/TIME.adt);
                        %QS(i)=min(QS(i),dnrev(i-1)*ds(i-1)*h0(i-1)/TIME.adt);        % <- basic approach, also works well, but bypass transport is a linear function of the time step
                    end
                end
                for i=[length(QS)-1:-1:1]
                    if QS(i)<0
                        QS(i)=max(QS(i),QS0(i+1)-dnrev(i)*ds(i)*h0(i)/TIME.adt);
                        %QS(i)=max(QS(i),-dnrev(i)*ds(i)*h0(i)/TIME.adt);        % <- basic approach, also works well, but bypass transport is a linear function of the time step
                    end
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
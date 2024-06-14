function [CHANNEL,COAST] = move_channel2(CHANNEL,COAST,WAVE,TIME)
% MOVE_CHANNEL when the width of the channels described by xr_mc, yr_mc
% is below channel_width the coast sections on either side are adapted by
% moving sand seaward according to the scheme dn -.5dn -.5dn
% where dn is fac*(actual width - channel_width) which is negative
% xr_mc    x coordinates of river(s) may be multiple separated by nan's
% yr_mc    y coordinates of river(s) may be multiple separated by nan's
% COAST.x_mc     x coordinates of coastline sections;
% y_mc     y coordinates of coastline sections;
% channel_width  vector with minimum width(s) of river(s)
% fac            factor (<1) applied to updating of coast sections
% x_inlet  x coordinate of most seaward point of channel where width is
%          less than channel_width
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

    ds0=COAST.ds0;
    if ~isscalar(ds0)
        ds0=min(ds0(:,3));
    end
    if CHANNEL.used~=0 && ~isempty(CHANNEL.xr_mc)
        plotdetails=0;
        ch='rgbm';
        section=get_sectionnr(COAST.x_mc); % for each point a section number is assigned
        %% Create grid along river axis
        
        [CHANNEL.xr_mc,CHANNEL.yr_mc]=initialize_grid(CHANNEL.xr_mc,CHANNEL.yr_mc,ds0);

        %% Center each river relative to the surrounding sections
        iriv=1;
        [ xr,yr,nr_mc,ir1(iriv),ir2(iriv) ] = get_one_polygon( CHANNEL.xr_mc,CHANNEL.yr_mc,iriv);
        [ x ,y ,n_mc ,i1(1),i2(1) ] = get_one_polygon( COAST.x_mc,COAST.y_mc,1);
        for iriv=1:nr_mc
            % Get one river
            [ xr,yr,nr_mc,ir1(iriv),ir2(iriv) ] = get_one_polygon( CHANNEL.xr_mc,CHANNEL.yr_mc,iriv);
            % Select the two closest distances to coast points
            nr=length(xr);
            if plotdetails
                figure(2)
                plot(COAST.x_mc,COAST.y_mc,'k',xr,yr,'b');
                hold on
            end
            for ir=1:nr
                if ir<nr
                    x1=xr(ir);x2=xr(ir+1);
                    y1=yr(ir);y2=yr(ir+1);
                else
                    x1=xr(ir-1);x2=xr(ir);
                    y1=yr(ir-1);y2=yr(ir);
                end
        %         d=(COAST.x_mc-x1)*(y2-y1)-(COAST.y_mc-y1)*(x2-x1);
        %         xleft=COAST.x_mc(d>0);yleft=COAST.y_mc(d>0);
        %         xright=COAST.x_mc(d<0);yright=COAST.y_mc(d<0);
                ext=200;
                xrpol=[xr(1)-ext*(yr(2)-yr(1)),xr,xr(end)-ext*(yr(end)-yr(end-1))];
                yrpol=[yr(1)+ext*(xr(2)-xr(1)),yr,yr(end)+ext*(xr(end)-xr(end-1))];
                inpol=inpolygon(COAST.x_mc,COAST.y_mc,xrpol,yrpol);
                xleft=COAST.x_mc(inpol);yleft=COAST.y_mc(inpol);
                xright=COAST.x_mc(~inpol);yright=COAST.y_mc(~inpol);
                [dmin1,imin1]=min(hypot(xr(ir)-xleft,yr(ir)-yleft));
                [dmin2,imin2]=min(hypot(xr(ir)-xright,yr(ir)-yright));

                x_d(ir,1)=xleft(imin1);
                y_d(ir,1)=yleft(imin1);
                x_d(ir,2)=xright(imin2);
                y_d(ir,2)=yright(imin2);
                ind1(ir)=find(COAST.x_mc==x_d(ir,1)&COAST.y_mc==y_d(ir,1),1);
                ind2(ir)=find(COAST.x_mc==x_d(ir,2)&COAST.y_mc==y_d(ir,2),1);
                if plotdetails
                    plot(xr(ir),yr(ir),'o',x_d(ir,1),y_d(ir,1),'r+',x_d(ir,2),y_d(ir,2),'g+')
                    hold on
                    plot(COAST.x_mc,COAST.y_mc,xrpol,yrpol)
                end
                %%pause
                alf=atan2d(y_d(ir,2)-y_d(ir,1),x_d(ir,2)-x_d(ir,1));
                b=hypot(y_d(ir,2)-y_d(ir,1),x_d(ir,2)-x_d(ir,1));
                a=0.5*(hypot(xr(ir)-x_d(ir,1),yr(ir)-y_d(ir,1))+ ...
                    hypot(xr(ir)-x_d(ir,2),yr(ir)-y_d(ir,2)));
                thetorg=atan2d(yr(ir)-y_d(ir,1),xr(ir)-x_d(ir,1))-alf;
                thet=acosd(.5*b/a)*sign(sind(thetorg));
                dxr=a*cosd(alf+thet);
                dyr=a*sind(alf+thet);
                % shift the river line to the middle
                xrnew(ir)=x_d(ir,1)+dxr;
                yrnew(ir)=y_d(ir,1)+dyr;
                %width(ir)=b;
                width(ir)=b/2;
            end
            [CHANNEL.xr_mc,CHANNEL.yr_mc]=insert_section(xrnew,yrnew,CHANNEL.xr_mc,CHANNEL.yr_mc,iriv);
            if plotdetails
                plot([xr',x_d(:,1)]',[yr',y_d(:,1)]','b')
                plot([xr',x_d(:,2)]',[yr',y_d(:,2)]','g')
                plot(xrnew,yrnew,[ch(iriv),'.'],'linewidth',1)
                axis equal
            end
            secnr1=find(diff(section(ind1))~=0);
            secnr2=find(diff(section(ind2))~=0);
            if isempty(secnr1)
                samedir1=ind1(end)>ind1(1)*ones(size(ind1));
            else
                samedir1(1:secnr1(1))=ind1(secnr1(1))>ind1(1);
                for i=1:length(secnr1)-1
                    samedir1(secnr1(i)+1:secnr1(i+1))=ind1(secnr1(i+1))>ind1(secnr1(i)+1);
                end
                samedir1(secnr1(end)+1:length(ind1))=ind1(end)>ind1(secnr1(end)+1);
            end
            if isempty(secnr2)
                samedir2=ind2(end)>ind2(1)*ones(size(ind2));
            else
                samedir2(1:secnr2(1))=ind2(secnr2(1))>ind2(1);
                for i=1:length(secnr2)-1
                    samedir2(secnr2(i)+1:secnr2(i+1))=ind2(secnr2(i+1))>ind2(secnr2(i)+1);
                end
                samedir2(secnr2(end)+1:length(ind2))=ind2(end)>ind2(secnr2(end)+1);
            end
            
            %% find inlet point (x_inlet,y_inlet)
            [minwidth,ind_inlet]=min(width);
            
            CHANNEL.x_inlet(iriv)=xrnew(ind_inlet);
            CHANNEL.y_inlet(iriv)=yrnew(ind_inlet);
            
            %% Find river discharge points on coast based on distance to inlet point
            dn=zeros(size(COAST.x_mc));
            if CHANNEL.disch_rate>0
                nmouthpts=3;
                dnmouth=CHANNEL.disch_rate/2/nmouthpts/ds0/COAST.h0*TIME.dt;
                ic1=ind1(ind_inlet);
                dn(ic1-1:ic1+1)=dnmouth;
                ic2=ind2(ind_inlet);
                dn(ic2-1:ic2+1)=dnmouth;
            end
            
            for ir=1:nr
                 if width(ir)<CHANNEL.width(iriv)/2
                    dw=CHANNEL.fac(iriv)*(width(ir)-CHANNEL.width(iriv)/2);
                    ic=ind1(ir);
                    [icm1,icm2,icp1,icp2]=find_neighbours(COAST.x_mc,COAST.y_mc,ic);
                    dn(ic)=dn(ic)+dw;
                    if samedir1(ir)
                        dn(icp1)=dn(icp1)-dw/2;
                        dn(icp2)=dn(icp2)-dw/2;
                    else
                        dn(icm1)=dn(icm1)-dw/2;
                        dn(icm2)=dn(icm2)-dw/2;
                    end
                    ic=ind2(ir);
                    [icm1,icm2,icp1,icp2]=find_neighbours(COAST.x_mc,COAST.y_mc,ic);
                    dn(ic)=dn(ic)+dw;
                    if samedir2(ir)
                        dn(icp1)=dn(icp1)-dw/2;
                        dn(icp2)=dn(icp2)-dw/2;
                    else
                        dn(icm1)=dn(icm1)-dw/2;
                        dn(icm2)=dn(icm2)-dw/2;
                    end
                end
            end
            for ic=1:length(COAST.x_mc)
                if dn(ic)~=0 
                    try
                    [icm1,icm2,icp1,icp2]=find_neighbours(COAST.x_mc,COAST.y_mc,ic);
                    catch
                        disp('problem in find_neighbours')
                    end
                    B=hypot(COAST.x_mc(icp1)-COAST.x_mc(icm1),COAST.y_mc(icp1)-COAST.y_mc(icm1));
                    dx(ic)=-dn(ic)*(COAST.y_mc(icp1)-COAST.y_mc(icm1))/B;
                    dy(ic)= dn(ic)*(COAST.x_mc(icp1)-COAST.x_mc(icm1))/B;
        %             COAST.x_mc(ic)=COAST.x_mc(ic)+dx(ic);
        %             COAST.y_mc(ic)=COAST.y_mc(ic)+dy(ic);
                end
            end
            for ic=1:length(COAST.x_mc)
                if dn(ic)~=0
                    COAST.x_mc(ic)=COAST.x_mc(ic)+dx(ic);
                    COAST.y_mc(ic)=COAST.y_mc(ic)+dy(ic);
                end
            end

            clear x_d; clear y_d;
            clear xrnew;clear yrnew;
        end
        %% Update grid along river axis
        if plotdetails
            plot(CHANNEL.xr_mc,CHANNEL.yr_mc,'.k','linewidth',1)
            plot(COAST.x_mc,COAST.y_mc,'m')
            %pause
        end

    end

end

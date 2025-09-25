function [COAST,TRANSP]=get_smoothangles(COAST,TRANSP)
% function [COAST,TRANSP]=get_smoothangles(COAST,TRANSP)
%
% limits the angles inbetween grid cells, 
% thus affecting spit width and stabilizing 
% small scale features in case of dense grids.
% 
% INPUT: 
%    COAST       : coastline data structure
%      .PHIc     : coastline orientation at each transport grid cell [°]
%      .maxangle : maximum coastline re-orientation between individual grid cells (default is 60°, with suggested range from 30° to 90°)
%    TRANSP      : transport data structure
%      .QS       : transport along the coast [m3/yr]
%
% OUTPUT:
%    COAST
%      .dPHIc    : computed angle difference between two subsequent coastline points
%    TRANSP
%      .QS       : transport along the coast [m3/yr]
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2021 IHE Delft & Deltares
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
    
    %% computed angle difference between two subsequent coastline points
    nw       = size(TRANSP.QS,1); % number of wave conditions taken along at this timestep (can be more than 1 in case of simultaneous wave conditions)
    nq       = size(TRANSP.QS,2); % number of coastline points
    COAST.dPHIc=mod(diff(COAST.PHIc)+180,360)-180;
    if COAST.clockwise && nw==1
        %% determine coastline points where angle is larger than the maxangle
        convexfactor=0.5; % factor determines sensitivity to max angles at concvex sections (e.g. if set to 0.5 then a gradual transition from 50% to 100% of maxangle is achieved for convex sections)
        idc1=find(COAST.dPHIc>COAST.maxangle*convexfactor & TRANSP.idrev(1:end-1)==0 & TRANSP.idrev(2:end)==0);
        concavefactor=0.5; % factor determines sensitivity to max angles at concave sections (e.g. if set to 0.5 then a gradual transition from 50% to 100% of maxangle is achieved for concave sections)
        idc2=find(COAST.dPHIc<-COAST.maxangle*concavefactor & TRANSP.idrev(1:end-1)==0 & TRANSP.idrev(2:end)==0);    
        idc=unique([idc1,idc2]);
        idc=setdiff(idc,find(TRANSP.shadow & abs(COAST.dPHIc)<90)); % exclude shadowed areas
        if ~COAST.cyclic
        idc=setdiff(idc,[1,2,length(COAST.dPHIc)-1,length(COAST.dPHIc)]); % exclude boundary points
        end
        
        if ~isempty(idc) 
            %% perform the smoothing of the transports
            QS0=TRANSP.QS;QS0a=QS0;
            % QS0s=get_smoothdata(QS0a,'');
            % QS0(QS0==0)=QS0s(QS0==0);
            QS1=TRANSP.QS;
            QS2=TRANSP.QS;
            factor = min(max(abs(2*abs(COAST.dPHIc(idc))-COAST.maxangle/2)/(COAST.maxangle/2),0),1);  % scale from 0 to 1 for the range +0.5*maxangle° to +1*maxangle°

            %% only perform smoothing when a disturbance is actually growing (either a spike or dent)
            ids=COAST.dPHIc(idc).*(QS0(idc)-QS0(idc+1))<0;
            factor(ids)=0;
            QSfac1=(QS0(idc)+QS0(idc+1))/2<=0; % smoothing only downdrift of the spike (to the left for negative transport)
            QSfac2=(QS0(idc)+QS0(idc+1))/2>=0; % smoothing only downdrift of the spike (to the right for positive transport)
            fac=max((QSfac1+QSfac2),1);
            QSfac1=QSfac1./fac;
            QSfac2=QSfac2./fac;

            %% smoothing at left side (transport point) of coastline point with large change in angle
            %QS1(idc)=(1-factor).*QS0(idc)+factor.*min(QS0(idc),QS0(idc+1)); % only transport in negative direction
            for nn=length(idc):-1:1 % sweep from right to left
                %QS1(idc(nn))=(1-factor(nn).*QSfac1(nn)).*QS0(idc(nn))+factor(nn).*QSfac1(nn).*QS1(idc(nn)+1);
                QS1(idc(nn))=(1-factor(nn).*QSfac1(nn)).*QS0(idc(nn))+factor(nn).*QSfac1(nn).*min(QS1(idc(nn)+1),QS0(idc(nn)));
            end
            
            %% smoothing at right side (transport point) of coastline point with large change in angle
            %QS2(idc+1)=(1-factor).*QS0(idc+1)+factor.*max(QS0(idc),QS0(idc+1)); % only transport in positive direction
            for nn=1:length(idc) % sweep from left to right
                %QS2(idc(nn)+1)=(1-factor(nn).*QSfac2(nn)).*QS0(idc(nn)+1)+factor(nn).*QSfac2(nn).*QS2(idc(nn));
                QS2(idc(nn)+1)=(1-factor(nn).*QSfac2(nn)).*QS0(idc(nn)+1)+factor(nn).*QSfac2(nn).*max(QS2(idc(nn)),QS0(idc(nn)+1));
            end
            
            %% combine smoothed transport at left and right side 
            ids=zeros(1,nq);
            ids(idc)=-factor.*QSfac1;
            ids(idc+1)=ids(idc+1)+factor.*QSfac2;
            ids=(ids+1)/2;
            TRANSP.QS=QS1.*(1-ids)+QS2.*ids;
            
            %% debug plot
            if 0
                idc1=idc((factor.*QSfac1)>0.3);
                idc2=idc((factor.*QSfac2)>0.3)+1;
                idc3=unique([idc1,idc2]);
                hold on;plot(COAST.xq(idc1),COAST.yq(idc1),'g*')
                hold on;plot(COAST.xq(idc2),COAST.yq(idc2),'gs')
                fprintf('%5.0f : %5.0f -> %5.0f \n',[idc3',QS0a(idc3)',TRANSP.QS(idc3)']');
                
                hold on;
                plot(COAST.x,COAST.y,'m-')
                plot(COAST.x(idc),COAST.y(idc),'m*')
                figure(111);clf;
                plot([1:length(COAST.x)+1]-0.5,QS0,'b-');hold on
                %plot([1:length(COAST.x)+1]-0.5,QS0s,'y-.');
                plot(idc-0.5,QS0(idc),'m*');
                plot(idc+0.5,QS0(idc+1),'m*');
                plot([1:length(COAST.x)+1]-0.5,TRANSP.QS,'r-')
                figure(11);
            end
        end
    end
end

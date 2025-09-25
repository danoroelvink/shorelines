function [TRANSP]=get_upwindcorrection(COAST,WAVE,TRANSP)
% function [TRANSP]=get_upwindcorrection(COAST,WAVE,TRANSP)
%
% The upwind correction function corrects the transport to the maximum transport 
% (as computed by the angles function with an S-Phi curve made for refracted nearshore breaking waves)
% at locations where the coastline orientation is larger than the critical high angle instability angle,
% with the requirement that the previous cell was still below the critical angle (with 'ref' suffix).
% 
% INPUT:  
%    COAST
%          .n           : Number of COAST.x points
%          .nq          : Number of TRANSP.QS points
%          .cyclic      : Index for COAST.cyclic coastline (0 or 1)
%    WAVE
%          .dPHItdp     : Relative angle of offshore waves with respect to the coast orientation ([Nx1] Radians)
%          .dPHIcrit    : Critical orientation of the coastline ( [Nx2] Radians)
%    TRANSP
%         option=0 : TRANSP.QS(1)=TRANSP.QSmax & Qs(2)=0; 
%         option=1 : Qs(1)=TRANSP.QSmax; Qs(2)=TRANSP.QSmax/2; Qs(3)=0; 
%         option=2 : Qs(1)=TRANSP.QSmax; Qs(2)=min(TRANSP.QS(2),TRANSP.QSmax/2); Qs(3)=0;
%          .twopoints   : Method for changes in sediment transport downdrift from the point with the upwind correction (i.e. at head of the spit) 
%          .shadowS     : Index of cells which are in the shadow zone due to coastline curvature (TRANSP.QS-points)
%          .shadowS_h   : Index of cells which are in the shadow zone due to hard structures (TRANSP.QS-points)
%          .QS          : Transport rates in grid cells [1xN] (in [m3/yr] including pores)
%          .QSmax       : Maximum transport for considered cells [1xN] (in [m3/yr] including pores)
%
% OUTPUT:
%    TRANSP
%          .QS          : Transport rates in grid cells [1xN] (in [m3/yr] including pores)
%          .im3         : Indices of coasltine points where an upwind correction was made with positive transport (for debugging)
%          .ip3         : Indices of coasltine points where an upwind correction was made with negative transport (for debugging)
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

    debuginfo=0; % use 1 or 2 for various output during debugging
    maxangle=max(COAST.maxangle,60);
    if isstruct(TRANSP.twopoints)
        TRANSP.twopoints=S.TRANSP.twopoints; % ensuring backward compatibility
    end
    
    xq=COAST.xq;
    yq=COAST.yq;
    nq=COAST.nq;
    dsq=COAST.dsq;
    nw=size(TRANSP.QS,1); % number of wave conditions taken along at this timestep (can be more than 1 in case of simultaneous wave conditions)

    TRANSP.ivals=zeros(nw,nq);
    QS=TRANSP.QS;
    QS0=TRANSP.QS;
    QS1=QS0;
    QS2=QS0;
    QSmax=TRANSP.QSmax;
    dPHIo=WAVE.dPHIo;
    dPHItdp=WAVE.dPHItdp;
    dPHIcrit=WAVE.dPHIcrit;

    if size(dPHIcrit,1)<nw
    dPHIcrit=repmat(dPHIcrit,[nw,1]);
    end
    if size(dPHIcrit,2)==1
    dPHIcrit=repmat(dPHIcrit,[1,nq]);
    end
    
    if COAST.clockwise==1 && TRANSP.suppresshighangle~=1
        for cw=[1,-1]  % positive and negative transport directions
            if COAST.cyclic
                iirange=[1:nq];
            else
                iirange=[2:nq-1];
            end
            
            % compute relative difference between 'coast angle' and 'critical coast angle'
            % both for the current TRANSP.QS-point ('i') and the previous one (left 'im1' or right 'ip1')
            dPHIcr1 = mod(dPHItdp-dPHIcrit+180,360)-180;          % index showing whether currrent point is high-angle in case of positive transport direction
            dPHIcr2 = mod(dPHItdp+dPHIcrit+180,360)-180;           % index showing whether currrent point is high-angle in negative transport direction
            dPHIcrSIGN = sign(dPHItdp).*(dPHItdp>dPHIcrit).*(abs(dPHItdp)-dPHIcrit) - sign(dPHItdp).*(dPHItdp<-dPHIcrit).*(abs(dPHItdp)-dPHIcrit); 
            
            dPHIcr=dPHIcr1;
            if cw==-1
                dPHIcr=dPHIcr2;
                % similar variable with just the locations where dPHItdp exceeds the dPHIcrit as non zero
                % in principle this 'dPHIcr' variable may replace the 4 dPHIcr's mentioned just above, but still needs to be tested and verified first
                iirange=fliplr(iirange);
            end
            
            for i=iirange
                if COAST.cyclic
                    im1=get_mod(i-1*cw,nq); % note that cw flips the indices for negative transport (cw=-1)
                    im2=get_mod(i-2*cw,nq); % note that cw flips the indices for negative transport (cw=-1)
                    ip1=get_mod(i+1*cw,nq);
                    ip2=get_mod(i+2*cw,nq);
                else     
                    if cw==1 % swap im1 and ip1 etc if direction of transport is negative (cw=-1)
                        im1=max(i-1,1);
                        im2=max(i-2,1);
                        ip1=min(i+1,nq);
                        ip2=min(i+2,nq);
                    elseif cw==-1 % swap im1 and ip1 etc if direction of transport is negative (cw=-1)
                        ip1=max(i-1,1);
                        ip2=max(i-2,1);
                        im1=min(i+1,nq);
                        im2=min(i+2,nq);
                    end
                end

                for kk=1:nw
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% Take into account the inertia of the flow, which means that transport does not decelerate/accelerate instanteneously but over a relaxation distance
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % the factor scales between 0 and 1 for accounting the difference in transport of the updrift cell with im1 index
                    % the relaxation length scales with the wave height
                    if ~isempty(TRANSP.relaxationlength)
                        relaxationlength=TRANSP.relaxationlength; %*min(WAVE.HStdp(im1),1);
                        relaxationmethod=2;
                        if relaxationmethod==1
                            % single (next downdrift) cell approach
                            fac=min(max(1-dsq(min(i,length(dsq)))/relaxationlength,0),1);
                            if QS0(kk,im1)>0 && QS0(kk,i)>=0 && cw==1
                                QS(kk,i)=max(QS(kk,i),QS0(kk,i)+max(QS0(kk,im1)-QS0(kk,i),0)*fac);
                            elseif QS0(kk,im1)<0 && QS0(kk,i)<=0 && cw==-1
                                QS(kk,i)=min(QS(kk,i),QS0(kk,i)+min(QS0(kk,im1)-QS0(kk,i),0)*fac);
                            end
                        elseif relaxationmethod==2
                            % multi-cell approach
                            fac0=dsq(min(i,length(dsq)))/relaxationlength;
                            fac1=[1:-fac0:0];
                            ind=min(max(i-cw*round([0:1:1/fac0]),1),nq);
                            if QS0(kk,im1)>0 && QS0(kk,i)>=0 && cw==1
                                QS1(kk,i)=max(QS0(kk,i)+max(QS1(kk,ind)-QS0(kk,i),0).*fac1);
                            elseif QS0(kk,im1)<0 && QS0(kk,i)<=0 && cw==-1
                                QS2(kk,i)=min(QS0(kk,i)+min(QS2(kk,ind)-QS0(kk,i),0).*fac1);
                            end
                            if cw==-1
                                QS(kk,i)=QS1(kk,i)+QS2(kk,i)-QS0(kk,i);
                            end
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% UPWIND CORRECTION FOR TRANSITION POINTS TO HIGH-ANGLE                                    %%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if cw*dPHIcr(kk,i)>0 && ~TRANSP.idrevq(i) ... % only when the wave angle is larger than the critical wave angle
                       && (cw*dPHIcr(kk,im1)<=0 || ...                                      % option 1: When there is a transition from low-angle to high-angle (so at 'im1' it is still low-angle)                        %(cw*QS(im1)>0 && max(cw*dPHIcr(kk,ip1),0)<max(cw*dPHIcr(kk,im1),0) && cw*dPHIcr(kk,i)>max(cw*dPHIcr(kk,im1),0)) || ...  % option 2: if it is just a local dip in the coastline, so the cell before it and after it are all dPHIcr>0
                          (cw*QS(kk,im1)>0 && (cw*dPHIcr(kk,im1)-cw*dPHIcr(kk,i))<-maxangle)) ... % option 3: When the updrift transport QS(im1) is larger than 0 and towards the cell & the angle change is very large
                       && ~TRANSP.shadowS(kk,im1) ...                                       % do not perform upwind for regions shaded by the coastline
                       && (isempty(TRANSP.shadowS_h)||(~isempty(TRANSP.shadowS_h)&&~TRANSP.shadowS_h(kk,ip1)&&~TRANSP.shadowS_h(kk,ip2))) ... % do not perform upwind for regions shaded by hard structures
                       && TRANSP.idrev(i)==0 && TRANSP.idrev(min(max(i-cw,1),nq))==0  % do not perform upwind for regions with revetments
                        
                        i1=i;                               % index of considered transition point (which is the first point that becomes 'high-angle')
                        i2=im1;                             % index of point updrift (which is still normal 'low-angle')                       
                        factorQS=1.05; % amplification factor of the max transport at the upwind corrected point, which smoothes local shapes if larger than 1
                        QS(kk,i1)=cw*max(factorQS*QSmax(kk,[i2,i1,ip1]));
                        cf=1.1;  % reduction factor of the transport downdrift of the upwind corrected point, which smoothes local shapes if smaller than 1
                        
                        if TRANSP.twopoints==1
                            QS(kk,ip1)=QS(kk,ip1)+cf*cw*max(0,(cw*QS(kk,i1)-cw*QS(kk,ip1))/2);                %cw*max(0,cw*QS(ip1));
                        elseif TRANSP.twopoints==2
                            QS(kk,ip1)=QS(kk,ip1)+cf*cw*max(0,(cw*QS(kk,i1)-cw*QS(kk,ip1))/2);                %cw*max(0.5*cw*QS(i1),cw*QS(ip1));
                            QS(kk,ip2)=QS(kk,ip2)+cf*cw*max(0,(cw*QS(kk,ip1)-cw*QS(kk,ip2))/2);               %cw*max(0,cw*QS(ip2));               
                        else
                            %QS(kk,ip1)=0;  
                        end
                        TRANSP.ivals(kk,i1)=cw; 
                        
                        % take into account the upwind correction for the smoothing of transport over the relaxation distance
                        if cw==1
                            QS1(kk,i1)=QS(kk,i1);
                            QS1(kk,ip1)=QS(kk,ip1);
                            QS1(kk,ip2)=QS(kk,ip2);
                        else
                            QS2(kk,i1)=QS(kk,i1);
                            QS2(kk,ip1)=QS(kk,ip1);
                            QS2(kk,ip2)=QS(kk,ip2);
                        end
                        
                        if debuginfo==1 
                            figure(100+cw);clf;
                            hold on;
                            plot(COAST.x,COAST.y,'k-');hold on;
                            hf=fill(COAST.x,COAST.y,'y');set(hf,'FaceColor',[1 1 0.5]);
                            plot(xq([i]),yq([i]),'ks');
                            plot(xq([i2]),yq([i2]),'ks');
                            plot(xq([i1]),yq([i1]),'ko');
                            plot(xq([ip1,ip2]),yq([ip1,ip2]),'k+');
                            xlim(xq([i])+[-1000,1000]);ylim(yq([i])+[-1000,1000]);
                            title('debug plot upwind scheme')
                        elseif debuginfo==2 
                            fprintf('-----------------------------------------------------------------\n');
                            ivalsorted=sort([i2,i1,ip1,ip2]);
                            fprintf('name       : %8s %8s %8s %8s  \n','i2','i1','ip1','ip2');
                            fprintf('index      : %8.0f %8.0f %8.0f %8.0f  \n',ivalsorted);
                            fprintf('QSmax      : %8.2f %8.2f %8.2f %8.2f  10^3 m^3/yr\n',QSmax(kk,ivalsorted)/10^3);
                            fprintf('QS         : %8.2f %8.2f %8.2f %8.2f  10^3 m^3/yr\n',QS(kk,ivalsorted)/10^3);
                            fprintf('dPHIo      : %8.3f %8.3f %8.3f %8.3f  °\n',dPHIo(kk,ivalsorted));
                            fprintf('dPHItdp    : %8.3f %8.3f %8.3f %8.3f  °\n',dPHItdp(kk,ivalsorted));
                            fprintf('dPHIo-tdp  : %8.3f %8.3f %8.3f %8.3f  °\n',dPHIo(kk,ivalsorted)-dPHItdp(kk,ivalsorted));
                            fprintf('dPHIcr>THR : %8.3f %8.3f %8.3f %8.3f  °\n',dPHIcr(kk,ivalsorted));
                            fprintf('QS(new)    : %8.2f %8.2f %8.2f %8.2f  10^3 m^3/yr\n',QS(kk,ivalsorted)/10^3);
                            %[QSmax(i2:ip2)/10^5;QS(kk,i2:ip2)/10^5;dPHIo(i2:ip2);dPHItdp(i2:ip2);dPHIo(i2:ip2)-dPHItdp(i2:ip2);dPHIcr(i2:ip2);i2:ip2]
                        end
                    end
                end
            end
        end
    end
    
    % debug value
    TRANSP.debug.QS2=QS;
    TRANSP.QS=QS;

end

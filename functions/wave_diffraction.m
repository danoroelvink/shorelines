function [STRUC,WAVE] = wave_diffraction(STRUC,COAST,WAVE,GROYNE,RUNUP)
% function [STRUC,WAVE] = wave_diffraction(STRUC,COAST,WAVE,GROYNE,RUNUP)
% 
% This routines computes the wave diffraction and transmission at groynes, offshore breakwaters and revetments. 
% The following approaches for diffraction can be used to compute the diffraction effect on the kd and direction.
% - kdform : diffraction approach for the wave height reduction (either 'Kamphuis' and 'Roelvink')
% - wdform : diffraction approach for the directional spreading of the waves (either 'Roelvink' of 'Dabees')
% The wave transmission is computed at offshore breakwaters. 
% 
% INPUT: 
%    STRUC
%        .diffraction        : switch for using diffraction (0/1)
%        .xhard              : x-coordinate of hard structures [m]
%        .yhard              : y-coordinate of hard structures [m]
%        .nrevet             : number of revetments (these are the last elements of xhard/yhard)
%        .kdform             : diffraction approach for the wave height reduction (either 'Kamphuis' or 'Roelvink')
%        .wdform             : diffraction approach for the directional spreading of the waves (either 'Roelvink' or 'Dabees')
%        .diffdist           : default maximum distance to wave points [m]
%        .wetstr_mc          : index indicating whether the structures are in the water or on land
%        .idgroyne           : index for each structure, indicating whether it is a groyne (with value = groyne number) or offshore breakwater/revetment (with value = 0) . 
%        .transmission       : switch for using wave transmission (0/1)
%        .transmbwdepth      : depth at breakwater location [m]
%        .transmcrestheight  : breakwater crest height [m]
%        .transmcrestwidth   : breakwater crest width [m]
%        .transmslope        : breakwater slope [-]
%        .transmform         : used formulation for calculating transmission ('angr', 'gent' or 'seabrhall')
%        .transmd50          : D50 of breakwater armour material [m]
%        .transmdir          : switch to determine wave propagation direction shoreward of breakwater (0/1)
%    COAST
%        .xq                 : x-coordinate of transport points (only current section)
%        .yq                 : y-coordinate of transport points (only current section)
%    WAVE
%        .PHItdp             : wave direction at the toe of the dynamic profile [°N]
%        .HStdp              : wave height at the toe of the dynamic profile [m]
%        .TP                 : wave period [s]
%        .dirspr             : directional spreading [°N]
%        .dist               : distance to the nearshore wave boundary points [m] (used as upper limit)
%    GROYNE
%        .phigroyne          : average orientation of the groyne (i.e. begin and end point compared to average of rest of groyne) [°N]
%    RUNUP
%        .swl                : still water-level due to combination of run-up, surge and tide.
% 
% OUTPUT:
%    STRUC
%        .xtip               : x-coordinates of diffraction points [m]
%        .ytip               : y-coordinates of diffraction points [m]
%        .Htip               : wave height at structure tips [m]
%        .phistruc           : orientation of structures [°N]
%        .xp                 : x-coordinate of transport points (only current section) [m]
%        .yp                 : y-coordinate of transport points (only current section) [m]
%    WAVE
%        .PHItdp             : wave direction after diffraction [°N]
%        .HStdp              : wave height after diffraction [m]
%        .dirtip             : wave directions at structure tips [°N]
%        .diff               : index with alongshore points with wave diffraction (summed up for all indices, for QS-points)
%        .diffxy             : index with alongshore points with wave diffraction (summed up for all indices, for xy-points)
%        .diff2              : index with alongshore points with wave diffraction (per coastal structure, for QS-points)
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2020 IHE Delft & Deltares
%
%       Dano Roelvink, Ahmed Elghandour
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
    
    STRUC.xq    = COAST.xq;    
    STRUC.yq    = COAST.yq;   
    WAVE.diff   = zeros(size(WAVE.PHItdp));  % wave diffraction at QS points
    % WAVE.diffxy = [];  % wave diffraction at XY points
    % WAVE.diff2  = {};    
    % WAVE.dirtip = [];
    % STRUC.xtip  = []; 
    % STRUC.ytip  = []; 
    % STRUC.Htip  = []; 
    epsom=45; % sheltering angle from sides of structures determining the size of the area where diffraction is computed.    

    if STRUC.diffraction==1 && ~isempty(STRUC.xhard)
        maxanglerotation=90.1;
        if ~isfield(WAVE,'dirspr')
            WAVE.dirspr=12;
        elseif isempty(WAVE.dirspr)
            WAVE.dirspr=12;
        end
        alfa=min(max(WAVE.dirspr-12,0),32)/20; 
        rotfac_longcrested = 0.8;   % at dirspr = 12° (default)
        omegat_longcrested = -20;    % at dirspr = 12° (default)
        rotfac_shortcrested = 0.8;  % at dirspr = 32° 
        omegat_shortcrested = -35;   % at dirspr = 32° 
        
        %% Get rotation factor and omegat (offset of angle at edge of influence area) of the diffraction routine
        wdform='Roelvink';
        rotfac = (1-alfa)*rotfac_longcrested + alfa*rotfac_shortcrested;
        omegat = (1-alfa)*omegat_longcrested + alfa*omegat_shortcrested;      
        if ~isempty(findstr(lower(STRUC.wdform),'dabees')) || ~isempty(findstr(lower(STRUC.wdform),'hurst')) || ~isempty(findstr(lower(STRUC.wdform),'kamph'))
            wdform='Dabees';
            omegat=0;
            rotfac=1;
        end
        
        %% Get input
        xhard=[];
        yhard=[];
        if ~isempty(STRUC.xhard) 
            xhard=[xhard,nan,STRUC.xhard];
            yhard=[yhard,nan,STRUC.yhard];
        end
        %if ~isempty(TRANSP.xsedlim) 
        %    xhard=[xhard,nan,TRANSP.xsedlim];
        %    yhard=[yhard,nan,TRANSP.ysedlim];
        %end
        if ~isempty(xhard)
            xhard=xhard(2:end);
            yhard=yhard(2:end);
        end
        
        x      = COAST.x;
        y      = COAST.y;
        xq     = COAST.xq;
        yq     = COAST.yq;
        nq     = length(xq);
        nw     = size(WAVE.PHItdp,1);
        
        for qq=1:nw
            PHItdp = WAVE.PHItdp(qq,:);
            HStdp  = WAVE.HStdp(qq,:);
            TP     = WAVE.TP(qq,:);  
            HStdp0 = HStdp;
            kdtdp  = ones(size(HStdp));
            
            %% Determine maximum distance from structure to diffraction points
            idnearwaves=[];
            if ~isempty(STRUC.diffdist)
                diffdist=STRUC.diffdist; % use predefined upper limit S.diffdist
            elseif ~isempty(WAVE.dist)
                diffdist=2*median(WAVE.dist); % use distance to wave locations as an upper limit
                %diffdist=COAST.s(end)/length(WAVE.WVC);
            else 
                diffdist=10000; % default upper limit of distance from structure to diffraction points
            end
            
            %% Determine diffraction points xtip,ytip for all structures
            [ x_struc,y_struc,nhard,~,~ ] = get_one_polygon( xhard,yhard,1);
            revetment=[];
            for j=1:nhard
                revetment(j)=0;
                [ x_struc,y_struc,nhard ] = get_one_polygon( xhard,yhard,j );
                [ wetstr ] = get_one_polygon( STRUC.wetstr_mc,j );
                
                % if required re-evaluate the wetstr
                % which may be necessary if an ofshore breakwater converts to a groyne
                if length(wetstr)~=length(x_struc)
                    nmc=length(COAST.x_mc)-1;
                    for ist=1:length(STRUC.xhard)
                        [~,icl]=min(hypot(COAST.x_mc-STRUC.xhard(ist),COAST.y_mc-STRUC.yhard(ist)));
                        if ~isnan(STRUC.xhard(ist))
                            im1=max(icl-1,1);
                            ip1=min(icl+1,nmc+1);
                            if isnan(COAST.x_mc(im1)); im1=im1-1;if im1==0;im1=2;end; end
                            if isnan(COAST.x_mc(ip1)); ip1=ip1+1;if ip1>nmc+1;ip1=nmc;end; end
                            dirm=360-atan2d(COAST.y_mc(ip1)-COAST.y_mc(im1),COAST.x_mc(ip1)-COAST.x_mc(im1)); 
                            dirstr=atan2d(STRUC.xhard(ist)-COAST.x_mc(icl),STRUC.yhard(ist)-COAST.y_mc(icl));
                            STRUC.wetstr_mc(ist)=int8(cosd(dirstr-dirm)>0);  % for octave, logicals cannot be assigned NaNs later on
                        else
                            STRUC.wetstr_mc(ist)=int8(0);
                        end
                    end
                    STRUC.wetstr_mc(isnan(STRUC.xhard))=nan;
                    [ wetstr ] = get_one_polygon( STRUC.wetstr_mc,j );
                end
                
                % Determine length along structure
                s_struc=[0,cumsum(hypot(diff(x_struc),diff(y_struc)))];
                mindist=0.2*s_struc(end);
                % Determine wave direction at the structure (one value based on
                % mean coordinates of structure)
                [~,closest]=min(hypot(mean(x_struc)-xq,mean(y_struc)-yq));
                PHItip=PHItdp(closest);
                % Project structure coordinates on line in wave direction
                proj=(x_struc-x_struc(1))*cosd(PHItip)-(y_struc-y_struc(1))*sind(PHItip);
                % The minimum and maximum distances to this line are the diffraction points
                
                %% determine left tip of structure
                proj(wetstr==0)=1e10;
                itip1=find(proj==min(proj));
                if length(itip1)>1
                    dist=[];
                    for dd=1:length(itip1)
                    dist(dd)=min(hypot(x_struc(itip1(dd))-xq,y_struc(itip1(dd))-yq));
                    end
                    itip1=itip1(dist==min(dist));
                    itip1=itip1(1);
                end
                xtip(1,j)=x_struc(itip1);
                ytip(1,j)=y_struc(itip1);
                
                % Locally throw out points too close to tip1 (only for offshore breakwater)
                if STRUC.idgroyne(j)==0
                    tooclose=hypot(x_struc-xtip(1,j),y_struc-ytip(1,j))<mindist;
                    tooclose(itip1)=0;
                    x_struc=x_struc(~tooclose);
                    y_struc=y_struc(~tooclose);
                    wetstr=wetstr(~tooclose);
                    proj=proj(~tooclose);
                end
                
                %% determine right tip of structure
                proj(wetstr==0)=-1e10;
                itip2=find(proj==max(proj));
                if length(itip2)>1
                    dist=[];
                    for dd=1:length(itip2)
                    dist(dd)=min(hypot(x_struc(itip2(dd))-xq,y_struc(itip2(dd))-yq));
                    end
                    itip2=itip2(dist==min(dist));
                    itip2=itip2(1);
                end
                xtip(2,j)=x_struc(itip2);
                ytip(2,j)=y_struc(itip2);
                
                %% in case of a revetment use the first and last point only
                if j>nhard-STRUC.nrevet
                    revetment(j)=1;
                    xtip(:,j)=[x_struc(1),x_struc(end)];
                    ytip(:,j)=[y_struc(1),y_struc(end)];
                end
                
                %% Locally trow out points too close to tip2;
                if STRUC.idgroyne(j)==0
                    tooclose=hypot(x_struc-xtip(2,j),y_struc-ytip(2,j))<mindist;
                    tooclose(xtip(1,j)==x_struc&ytip(1,j)==y_struc)=0;
                    tooclose(itip2)=0;
                    x_struc=x_struc(~tooclose);
                    y_struc=y_struc(~tooclose);
                end
                % insert x_struc, y_struct back into xhard, yhard
                [xhard,yhard]=insert_section(x_struc,y_struc,xhard,yhard,j);
                % Get wave direction and wave height at each tip
                for tip=1:2
                    [~,closest]=min(hypot(xtip(tip,j)-xq,ytip(tip,j)-yq));
                    Htip(tip,j)=HStdp(closest); % used for wave diffraction / sheltering
                    TPtip(tip,j)=TP(closest); % used for wave diffraction / sheltering
                    dirtip(tip,j)=PHItdp(closest);
                    itip(tip,j)=closest(1);
                end
                Htip(3,j)=((Htip(1,j).^2+Htip(2,j).^2)/2).^0.5; % for wave transmission
                TPtip(3,j)=((TPtip(1,j).^2+TPtip(2,j).^2)/2).^0.5; % for wave transmission
                sinkd=(Htip(1,j).^2.*sind(dirtip(1,j)) + Htip(2,j).^2.*sind(dirtip(2,j)));
                coskd=(Htip(1,j).^2.*cosd(dirtip(1,j)) + Htip(2,j).^2.*cosd(dirtip(2,j)));
                dirtip(3,j)=mod(atan2d(sinkd,coskd),360); % direction affects the 'effective width' of the structure for wave transmission
                
                if 0
                    figure(88);clf;
                    %plot(COAST.x,COAST.y,'y');hold on;
                    plot(COAST.x_mc,COAST.y_mc,'y','LineWidth',4);hold on;
                    plot(xhard,yhard,'k-','LineWidth',0.5,'Color',[0.7 0.7 0.7]);
                    plot(x_struc,y_struc,'k-','LineWidth',2);
                    plot(xtip(1,:),ytip(1,:),'gs');
                    plot(xtip(2,:),ytip(2,:),'ms');
                    plot(x_struc(1),y_struc(1),'ks')
                end
                if 0
                    hold on;
                    plot(xtip(1,:),ytip(1,:),'gs');
                    plot(xtip(2,:),ytip(2,:),'ms');
                    plot(x_struc,y_struc,'k-','LineWidth',2);
                    plot(x_struc(1),y_struc(1),'ks')
                end
                
                %% Remove hard structures which are too far away
                xpp = [xtip(1,j):(xtip(2,j)-xtip(1,j))/9:xtip(2,j)];
                ypp = [ytip(1,j):(ytip(2,j)-ytip(1,j))/9:ytip(2,j)];
                if isempty(xpp);xpp=repmat(xtip(1,j),size(ypp));end
                if isempty(ypp);ypp=repmat(ytip(1,j),size(xpp));end
                idnearwaves(j,1:length(xq))=0;
                for pp=1:length(xpp)
                    distxq=((xq-xpp(pp)).^2+(yq-ypp(pp)).^2).^0.5;
                    distxq=get_smoothdata(distxq,'mean',10);
                    Lgroynes=(diff(GROYNE.x).^2+diff(GROYNE.y).^2).^0.5;
                    idpp=distxq<max(diffdist,mean(Lgroynes(:))*2);
                    if ~isempty(idpp)
                    idnearwaves(j,idpp)=1;
                    end
                end
                
                if STRUC.idgroyne(j)==0
                    STRUC.phistruc(j,1)=mod(atan2d(diff(xtip(:,j)),diff(ytip(:,j))),360);
                    STRUC.phistruc(j,2)=mod(atan2d(-diff(xtip(:,j)),-diff(ytip(:,j))),360);
                else
                    gg=STRUC.idgroyne(j);
                    STRUC.phistruc(j,1)=mod(GROYNE.phigroyne(gg),360);
                    STRUC.phistruc(j,2)=mod(GROYNE.phigroyne(gg),360);
                end
            end
            
            %% Compute angle to diffraction point for each transport point,
            %% relative to incident wave angle
            % phitip1 and phitip2 are angles from coast transport point to diffraction
            % points [nhard by nq]
            phitip1=mod(atan2d(xtip(1,:)'-xq,ytip(1,:)'-yq),360); % matrix [nhard by nq]
            phitip2=mod(atan2d(xtip(2,:)'-xq,ytip(2,:)'-yq),360);
            
            % omega1 and omega2 angles within diffraction zone, positive inwards
            %approach='new';
            omega1=[];
            omega2=[];
            for j=1:nhard
                if STRUC.idgroyne(j)==0
                    % original method
                    omega1(j,:)=repmat(dirtip(1,j)',1,nq)-phitip1(j,:);     % matrix [nhard by nq]
                    omega2(j,:)=phitip2(j,:)-repmat(dirtip(2,j)',1,nq);
                    omega1(j,:)=mod(omega1(j,:)+180,360)-180;
                    omega2(j,:)=mod(omega2(j,:)+180,360)-180;
                else
                    % benefit of this method is that later cross-sections with structures do not need to be computed anymore
                    phis1=repmat(STRUC.phistruc(j,1),[1,nq]);
                    phis2=repmat(STRUC.phistruc(j,2),[1,nq]);
                    PHIWAVEwrtGRO1=mod(repmat(dirtip(1,j)'-STRUC.phistruc(j,1),[1,nq])+180,360)-180;
                    PHIWAVEwrtGRO2=mod(repmat(STRUC.phistruc(j,2)-dirtip(2,j)',[1,nq])+180,360)-180;
                    PHICOASTwrtGRO1=mod(phitip1(j,:)-phis1+180,360)-180;
                    PHICOASTwrtGRO2=mod(180+phis2-phitip2(j,:),360)-180;
                    omega1(j,:)=PHIWAVEwrtGRO1-PHICOASTwrtGRO1;
                    omega2(j,:)=PHIWAVEwrtGRO2-PHICOASTwrtGRO2;
                end
            end
            
            % % omega1 and omega2 angles within diffraction zone, positive inwards
            % approach='new';
            % if strcmpi(approach,'old')
            %     % original method
            %     omega1b=repmat(dirtip(1,:)',1,nq)-phitip1;     % matrix [nhard by nq]
            %     omega2b=phitip2-repmat(dirtip(2,:)',1,nq);
            %     omega1b=mod(omega1b+180,360)-180;
            %     omega2b=mod(omega2b+180,360)-180;
            % elseif strcmpi(approach,'new')
            %     % benefit of this method is that later cross-sections with structures do not need to be computed anymore
            %     phis1=repmat(STRUC.phistruc(:,1),[1,nq]);
            %     phis2=repmat(STRUC.phistruc(:,2),[1,nq]);
            %     PHIWAVEwrtGRO1=mod(repmat(dirtip(1,:)'-STRUC.phistruc(:,1),[1,nq])+180,360)-180;
            %     PHIWAVEwrtGRO2=mod(repmat(STRUC.phistruc(:,2)-dirtip(2,:)',[1,nq])+180,360)-180;
            %     PHICOASTwrtGRO1=mod(phitip1-phis1+180,360)-180;
            %     PHICOASTwrtGRO2=mod(180+phis2-phitip2,360)-180;
            %     omega1b=PHIWAVEwrtGRO1-PHICOASTwrtGRO1;
            %     omega2b=PHIWAVEwrtGRO2-PHICOASTwrtGRO2;
            % end
            
            %% om is the change of direction due to diffraction
            om=zeros(size(HStdp));
            % om1, om2 are wave directions relative to incident direction due to each
            % diffraction point, rotated by % rotfac (=2 according to Hurst et al,
            % but 0.8 according to comparison with % wave-resolving XBeach model.
            % Both om1 and om2 are zero at edge of influence zone.
            om1=min(max(0,(omega1-omegat)*rotfac),maxanglerotation);
            om2=-min(max(0,(omega2-omegat)*rotfac),maxanglerotation);
            om3=zeros(size(om1));  % this assumes the wave direction does not change over the breakwater after wave transmission
            
            %% Diffraction coefficients according to simple formula
            % kdform = 'Roelvink' or 'Kamphuis'
            kdform = 'Roelvink';
            if ~isempty(findstr(lower(STRUC.kdform),'kamp'))
                kdform='Kamphuis';
            end
            kd1=wave_diffraction_coeff(omega1,kdform,wdform,WAVE.dirspr);
            kd2=wave_diffraction_coeff(omega2,kdform,wdform,WAVE.dirspr);
            
            %% Transmission over breakwater if submerged
            kd3=kd1.*0;
            for j=1:nhard           
                if STRUC.transmission==1 && STRUC.idgroyne(j)==0 %if wave propagation can take place over a submerged breakwater, then calculate amount of wave energy transmission
                    xhard=xtip(1:2,j);
                    yhard=ytip(1:2,j);
                    % calculate wave transmission
                    [kd3obw]=wave_transmission(xhard,yhard,xq,yq,x,y,RUNUP.swl,Htip(3,j),TPtip(3,j),STRUC.transmbwdepth(j),STRUC.transmcrestheight(j),STRUC.transmcrestwidth(j),STRUC.transmslope(j),STRUC.transmform{j},STRUC.transmd50(j)); %calculate amount of wave energy transmission 
                    kd3(j,:)=kd3obw;
                    % calculate wave propagation direction shoreward of breakwater - om3 
                    if STRUC.transmdir==1 % determines wave propagation direction shoreward of breakwater
                        PHIWAVEwrtGRO3=mod((dirtip(3,j)'-STRUC.phistruc(j,1)+90)+180,360)-180;
                        PHIWAVEwrtGRO3(2)=mod((dirtip(3,j)'-STRUC.phistruc(j,1)-90)+180,360)-180;
                        [~,ii]=min(abs(PHIWAVEwrtGRO3)); % use angle between incoming wave angle and breakwater
                        id3=find(abs(om1(j,:))>0 & abs(om2(j,:))>0);    
                        om3(j,id3)=PHIWAVEwrtGRO3(ii(1));
                        
                        % the rotation due to the wave transmission cannot be more than for the nearest diffraction point
                        om3(j,:)=sign(om3(j,:)).*min(abs(om3(j,:)),min(abs(om1(j,:)),abs(om2(j,:)))).*(1-kd3(j,:));
                        % alternative way of smoothing the om3 rotation due to transmission at the edges
                        % nsteps=max(200/COAST.ds0,1);
                        % om3(j,:)=get_smoothdata(om3(j,:),'',nsteps).*(1-kd3(j,:));
                    end 
                end
            end      
            
            %% Include wave transmission in the diffraction 
            % scales the kd1 and kd2 with respect to the kd3 coefficient
            % modifies the energy kd1/kd2 and direction of the om1/om2
            % treated per individual structure
            for j=1:nhard
                if STRUC.idgroyne(j)==0 
                    kd1_ratio = kd1(j,:).^2./(kd1(j,:).^2+kd2(j,:).^2);
                    kd2_ratio = kd2(j,:).^2./(kd1(j,:).^2+kd2(j,:).^2);
                    kd1(j,:)=(kd1(j,:).^2.*kd1_ratio.*(1-kd3(j,:).^2)).^0.5;
                    kd2(j,:)=(kd2(j,:).^2.*kd2_ratio.*(1-kd3(j,:).^2)).^0.5;  
                    
                    if sum(kd3(j,:))~=0
                        % include the wave transmission (kd3) in the om1/om2 and kd1/kd2 components
                        sinkd1=kd1(j,:).^2.*sind(om1(j,:))+kd1_ratio.*kd3(j,:).^2.*sind(om3(j,:));
                        coskd1=kd1(j,:).^2.*cosd(om1(j,:))+kd1_ratio.*kd3(j,:).^2.*cosd(om3(j,:));
                        om1(j,:)=atan2d(sinkd1,coskd1);
                        kd1(j,:)=(kd1(j,:).^2+kd1_ratio.*kd3(j,:).^2).^0.5;
                        %kd1(j,:)=(sinkd1+coskd1).^(1/2);
                        sinkd2=kd2(j,:).^2.*sind(om2(j,:))+kd2_ratio.*kd3(j,:).^2.*sind(om3(j,:));
                        coskd2=kd2(j,:).^2.*cosd(om2(j,:))+kd2_ratio.*kd3(j,:).^2.*cosd(om3(j,:));
                        om2(j,:)=atan2d(sinkd2,coskd2);
                        kd2(j,:)=(kd2(j,:).^2+kd2_ratio.*kd3(j,:).^2).^0.5;
                        %kd2(j,:)=(sinkd2+coskd2).^(1/2);
                        kd3(j,:)=0; % not further use kd3
                    end
                end
            end
            
            %% Determine which points are in which influence zone
            in=nan(nhard,size(COAST.xq,2));
            for j=1:nhard
                if STRUC.idgroyne(j)~=0 
                    % groyne
                    gg=STRUC.idgroyne(j);
                    
                    % determine from which side the wave comes
                    % do not use if wave comes from sheltered side of groyne (i.e. when GROYNEside(j)==1)
                    GROYNEside=sign(mod(dirtip(1,:)'-STRUC.phistruc+180,360)-180);     
                       
                    if GROYNE.idcoast(gg,1)>=COAST.i_mc %&& GROYNEside(j)~=1 
                        % coast left of groyne
                        in(j,omega1(j,:)>=(omegat-epsom))=-gg; 
                        kd2(j,:)=0;
                        kd3(j,:)=0;
                    elseif GROYNE.idcoast(gg,2)<=COAST.i_mc %&& GROYNEside(j)~=1 
                        % coast right of groyne
                        in(j,omega2(j,:)>=(omegat-epsom))=gg; 
                        kd1(j,:)=0;
                        kd3(j,:)=0;
                    end
                elseif revetment(j)==1
                    % revetment
                    % it is treated similar to the offshore breakwater, but the region of influence ('sahadow angle') is much smaller. 
                    in(j,omega1(j,:)>=(omegat-10)&omega2(j,:)>=(omegat-10))=0; % matrix [nhard by nq]
                    in(j,itip(1,j)+1:itip(2,j)-1)=nan; % exclude the revetment itself
                else
                    % offshore breakwater
                    %in(j,omega1(j,:)>=(omegat-epsom)&omega2(j,:)>=(omegat-epsom))=0; % matrix [nhard by nq]
                    in(j,:)=0;
                end
                % remove diffraction points that are too far from structure
                in(j,~idnearwaves(j,:))=nan;
            end
            
            %% Shadow the tips of the structures that are hidden behind others
            % This is elevant for constructions far behind others. 
            % For example this groyne behind a breakwater
            % The effect of the other structure on the kd1, kd2 and kd3 is included. 
            % Effects of shielding structure on wave rotation at the tip is not included yet!
            %
            %         xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            %
            %
            %                         X
            %                         X
            % ========================X=============================
            for jj=1:nhard
                % find whether the two tips and the middle of the breakwater cross with any structure
                xs=STRUC.xhard;
                ys=STRUC.yhard;
                xp=[xtip(:,jj);mean(xtip(:,jj))];
                yp=[ytip(:,jj);mean(ytip(:,jj))];
                [shadowS,jstruc]=find_shadows_mc(xp(:)',yp(:)',xs(:)',ys(:)',mod(dirtip(1:3,jj)',360));
                % exclude a crossing with itself
                if ~isempty(jstruc)
                for i=1:3
                    jstruc{i}=setdiff(jstruc{i},jj);
                    shadowS(i)=~isempty(jstruc{i});
                end
                end
                
                % shield left side from structures in jstruc{1}
                if shadowS(1)~=0
                    for cc=1:length(jstruc{1})
                        kd1(jj,:)=kd1(jj,:) .* (kd1(jstruc{1}(cc),:).^2+kd2(jstruc{1}(cc),:).^2+kd3(jstruc{1}(cc),:).^2).^0.5;
                        % % find area that is shielded by jstruc
                        % id1=(om2(jstruc{1},:)<0); % id1 is where the structure on the left side is the constraint for the directions
                        % id2=(om1(jstruc{1},:)>0);  % id2 is where the structure on the right side is the constraint for the directions
                        % id3=id1&id2;        % id3 is where they both are influenced, here take the average
                        % id1=id1&~id3;
                        % id2=id2&~id3;  
                        % % make sure that wave angles are limited by both sides of the gap
                        % om1(jj,id1)=min(om1(jj,id1),om2(jj2,id1));
                        % om2(jj,id2)=max(om1(jj,id2),om2(jj2,id2));
                        % om1(jj,id3)=(om1(jj,id3)+om1(jj2,id3))/2; 
                        % om2(jj,id3)=(om2(jj,id3)+om2(jj2,id3))/2; 
                    end
                end
                % shield right side from structures in jstruc{2}
                if shadowS(2)~=0
                    for cc=1:length(jstruc{2})
                    kd2(jj,:)=kd2(jj,:) .* (kd1(jstruc{2}(cc),:).^2+kd2(jstruc{2}(cc),:).^2+kd3(jstruc{2}(cc),:).^2).^0.5;
                    end
                end
                % shield middle section from structures in jstruc{3}
                if shadowS(3)~=0
                    for cc=1:length(jstruc{3})
                    kd3(jj,:)=kd3(jj,:) .* (kd1(jstruc{3}(cc),:).^2+kd2(jstruc{3}(cc),:).^2+kd3(jstruc{3}(cc),:).^2).^0.5;
                    end
                end           
            end       
            
            %% Compute relative impact of hard tips on each other, so how a tip is shaded by another tip (not used)
            if 0
                for jj=1:nhard
                    phitip1hard(jj,:)=mod(atan2d(xtip(1,jj)'-xtip(1,:),ytip(1,jj)'-ytip(1,:)),360); % matrix [nhard by nq]
                    phitip2hard(jj,:)=mod(atan2d(xtip(2,jj)'-xtip(2,:),ytip(2,jj)'-ytip(2,:)),360); % matrix [nhard by nq]
                end
                omega1tip=repmat(dirtip(1,:)',1,nhard)-phitip1hard;     % matrix [nhard by nq]
                omega2tip=phitip2hard-repmat(dirtip(2,:)',1,nhard);
                omega1tip=mod(omega1tip+180,360)-180;
                omega2tip=mod(omega2tip+180,360)-180;
                om1tip=min(max(0,(omega1tip-omegat)*rotfac),90);
                om2tip=-min(max(0,(omega2tip-omegat)*rotfac),90);
                kd1tip=wave_diffraction_coeff(omega1tip,kdform);
                kd2tip=wave_diffraction_coeff(omega2tip,kdform);
                for jj=1:size(xtip,2) 
                    kd1tip(jj,1:jj)=1;
                    kd2tip(jj,jj:end)=1;
                end
            end
            
            %% Discard the effects of structures that are outside the model area by only taking the nearest structure
            % Project the alongshore position of diffraction points of structures on the xq grid
            idstruc=nan(2,length(HStdp));
            idS=[];distS=[];typeS=[];
            for jj=1:nhard
                [~,~,idgrid,distw]=get_interpolation_on_grid('weighted_distance',xtip(:,jj),ytip(:,jj),xq,yq);
                idgrid=round(idgrid);
                idstruc(1,min(idgrid))=jj;
                idstruc(2,max(idgrid))=jj;
                idS(jj,:)=idgrid(:)';
                distS(jj,:)=distw(:)';
                if STRUC.idgroyne(jj)~=0
                    typeS(jj,:)=[0,0];
                    idstruc(1,min(idgrid))=-jj;
                    idstruc(2,max(idgrid))=-jj;
                else
                    typeS(jj,:)=[1e6,1e6];
                end
            end
            
            % The nearest groyne is prefered at each boundary. 
            % if no groynes are available, then the nearest offshore breakwater is chosen. 
            id1=find(idS(:,2)==1);
            id2=find(idS(:,1)==length(xq));
            if length(id1)>1 % left boundary
                dist2=distS+1e10;
                dist2(id1)=(distS(id1,1)+typeS(id1,1));
                idgroyne=find(dist2==min(dist2),1); % find nearest structure, often a groyne.
                idnot=setdiff(id1,idgroyne); % make sure to use only the structures at the left boundary (i.e. id1)   
                idS(idnot,:)=nan; 
                kd1(idnot,:)=0;
                kd2(idnot,:)=0;
                in(idnot,:)=nan;
            end
            if length(id2)>1 % right boundary
                dist2=distS+1e10;
                dist2(id2)=(distS(id2,1)+typeS(id2,1));
                idgroyne=find(dist2==min(dist2),1); % find nearest structure, often a groyne.
                idnot=setdiff(id2,idgroyne); % make sure to use only the structures at the right boundary (i.e. id2)
                idS(idnot,:)=nan; 
                kd1(idnot,:)=0;
                kd2(idnot,:)=0;
                in(idnot,:)=nan;
            end
            
            %% Correct kd values in case there is a gap between two structures with influence of both structures (either 2 breakwaters, breakwater with groyne or 2 groynes)
            % relevant for gaps between structures
            % for example the gap between this groyne and breakwater
            %
            %   xxxxxxxxxxxxxxxxxxxxxx             X
            %                                      X
            %                                      X
            %                                      X
            %                                      X
            % =====================================X================
            kd1plot=kd1;
            kd2plot=kd2;
            kd3plot=kd3;
            om1plot=om1;
            om2plot=om2;
            om3plot=om3;
            in(om2==-maxanglerotation)=nan; % 
            in(om1==maxanglerotation)=nan;  % 
            for jj=1:nhard
                notjj=setdiff(find(idS(:,2)~=0),jj);
                %% correct the diffraction in gaps between structures 
                % (kd1 of right structure and kd2 of the left structure)
                indexdiff=-(idS(notjj,2)-idS(jj,1));
                indexdiff(indexdiff<0)=nan; % only use structures to the left of the considered structure
                id=notjj(find(indexdiff==abs(min(indexdiff)),1));
                if ~isempty(id)
                    jj0=[jj;id];
                    in0=find(sum(~isnan(in(jj0,:)),1)>1); % only when the influence areas of at least two structures overlap!
                    in1=[idS(id,1):idS(jj,1)];
                    jj2=id(1);
                    
                    % make sure that wave energy through the gap is limited by the multiplication of both kd components (i.e. reduction by both)
                    kd1ratio=kd2(jj2,in0);
                    kd2ratio=kd1(jj,in0);
                    kd1(jj,in0)=kd1(jj,in0) .* kd1ratio /sqrt(2);     % correct left side kd1 
                    kd2(jj2,in0)=kd2(jj2,in0) .* kd2ratio /sqrt(2);   % correct right side kd2
                    
                    % make sure that wave angles are limited by both sides of the gap
                    %figure;plot(om1(jj,:),'b');hold on;plot(om2(jj2,:),'r');
                    id1=in0(om2(jj2,in0)<0); % id1 is where the structure on the left side is the constraint for the directions
                    id2=in0(om1(jj,in0)>0);  % id2 is where the structure on the right side is the constraint for the directions
                    id3=intersect(id1,id2);        % id3 is where they both are influenced, here take the average
                    id1=setdiff(id1,id3);
                    id2=setdiff(id2,id3);
                    
                    om1(jj,id1)=min(om1(jj,id1),om2(jj2,id1));
                    om2(jj2,id2)=max(om1(jj,id2),om2(jj2,id2));
                    om1(jj,id3)=(om1(jj,id3)+om2(jj2,id3))/2; %(om1(jj,id3)+om2(jj2,id3))/2;
                    om2(jj2,id3)=(om1(jj,id3)+om2(jj2,id3))/2; %om1(jj,id3);
                    %plot(om1(jj,:),'g:','LineWidth',2);plot(om2(jj2,:),'m--','LineWidth',1);
                end
            end
                    
            %% Use only the most dominantly impacting groyne or breakwater per side 
            % this routine is only used for groynes, where one groyne partly extents over another groyne that is closer to the coast. 
            % in this case part of the coastal section should be influenced just by groyne 1 (where dominant) and the other part of the beach just by groyne 2 (where that one is dominant).
            % So: kd1 and kd2 should have just one representative diffraction point and not two (or more). 
            % The diffraction point with the largest wave angle rotation is chosen. 
            % Whether it is also relevant for offshore breakwater combinations (with groynes) should be considered.
            for i=1:size(in,2)
                in2=find(~isnan(in(:,i)));
                % for groynes:
                in2neg=in2(in(in2,i)<0);
                in2pos=in2(in(in2,i)>0);
                id0=find(kd1(in2neg,i)>min(kd1(in2neg,i)));
                id1=find(kd2(in2pos,i)>min(kd2(in2pos,i)));
                if ~isempty(id0)
                in(in2neg(id0),i)=nan;
                end
                if ~isempty(id1)
                in(in2pos(id1),i)=nan;
                end
                % for offshore breakwaters:
                % in2is=in2(in(in2,i)==0);
                % id0=find(kd1(in2is,i)>min(kd1(in2is,i)),1);
                % id1=find(kd2(in2is,i)>min(kd2(in2is,i)),1);
                % in(id0,i)=nan;
                % in(id1,i)=nan;
            end
            
            %% DEBUG PLOTS
            if 0
                figure;
                plot(COAST.x_mc,COAST.y_mc);hold on;plot(STRUC.xhard,STRUC.yhard,'r');
                idnan=find(isnan(STRUC.xhard));
                if ~isempty(idnan)
                    idnan=[0,idnan,length(STRUC.xhard)+1];
                end
                for rr=1:length(idnan)-1
                    ids=[idnan(rr)+1:idnan(rr+1)-1];
                    text(mean(STRUC.xhard(ids)),mean(STRUC.yhard(ids)),num2str(rr));
                end
                figure;
                subplot(2,1,1);plot(WAVE.HStdp);
                subplot(2,1,2);plot(WAVE.PHItdp);
            end
                   
            %% check shielding of the coastline on wave diffraction points
            for j=1:nhard
                inn=~isnan(in(j,:)); % 1d mask for influence of structure j
                
                %% Check if diffraction points can be 'seen' from coast transport points
                % If not, setCOAST. kd to 0. epsx and epsy are to avoid counting
                % intersections in offshore parts of structure with finite width
                irange=find(inn);
                inshadow1=zeros(1,size(inn,2));
                inshadow2=zeros(1,size(inn,2));
                sd1=sind(PHItdp-om1(j,:));
                cd1=cosd(PHItdp-om1(j,:));
                sd2=sind(PHItdp-om2(j,:));
                cd2=cosd(PHItdp-om2(j,:));
                for i=irange
                    % only evaluate when there is a non-zero diffraction
                    if kd1(j,i)~=0 
                        % search for crossings with the coastline at the 'left' side (i.e. for 'kd1')
                        dx=abs(xq(i)-xtip(1,j));
                        dy=abs(yq(i)-ytip(1,j));
                        Ldiffpoint=(dx.^2+dy.^2).^0.5;
                        epsxc=0.001*sd1(i)*Ldiffpoint;
                        epsyc=0.001*cd1(i)*Ldiffpoint;
                        [xx_coast1,yy_coast1]=get_intersections([xq(i)+epsxc xq(i)+sd1(i)*Ldiffpoint*0.8],[yq(i)+epsyc yq(i)+cd1(i)*Ldiffpoint*0.8],xq(irange),yq(irange));   
                        %epsxc=0.001*(xtip(1,j)-xq(i));
                        %epsyc=0.001*(ytip(1,j)-yq(i));
                        %[xx_coast1,yy_coast1]=get_intersections([xq(i)+epsxc xtip(1,j)-epsxc/1000],[yq(i)+epsyc ytip(1,j)-epsyc/1000],xq(irange),yq(irange));   
                        if ~isempty(xx_coast1)
                            inshadow1(i)=1;
                        end
                    end
                    if kd2(j,i)~=0
                        % search for crossings with the coastline at the 'right' side (i.e. for 'kd2')
                        dx=abs(xq(i)-xtip(2,j));
                        dy=abs(yq(i)-ytip(2,j));
                        Ldiffpoint=(dx.^2+dy.^2).^0.5;
                        epsxc=0.001*sd2(i)*Ldiffpoint;
                        epsyc=0.001*cd2(i)*Ldiffpoint;
                        [xx_coast2,yy_coast2]=get_intersections([xq(i)+epsxc xq(i)+sd2(i)*Ldiffpoint*0.8],[yq(i)+epsyc yq(i)+cd2(i)*Ldiffpoint*0.8],xq(irange),yq(irange));   
                        %epsxc=0.001*(xtip(2,j)-xq(i));
                        %epsyc=0.001*(ytip(2,j)-yq(i));
                        %[xx_coast2,yy_coast2]=get_intersections([xq(i)+epsxc xtip(2,j)-epsxc/1000],[yq(i)+epsyc ytip(2,j)-epsyc/1000],xq(irange),yq(irange));
                        if ~isempty(xx_coast2)
                            inshadow2(i)=1;
                        end
                    end
                end
    
                if revetment(j)==0
                    % use shadowing of coastline if it is more than 1 cell
                    % and use 50% reduction at the first shadowed coastal cell and 100% reduction in subsequent cells
                    % for the left side of the structure 'kd1'
                    b=inshadow1(1:end-1)==1 & inshadow1(2:end)==1;
                    c=[(b(1:end-1)+b(2:end))/2];
                    inshadow1i=[c(1),c,c(end)];
                    kd1(j,:)=kd1(j,:).*(1-inshadow1i);
    
                    % use shadowing of coastline if it is more than 1 cell
                    % and use 50% reduction at the first shadowed coastal cell and 100% reduction in subsequent cells
                    % for the right side of the structure 'kd2'
                    b=inshadow2(1:end-1)==1 & inshadow2(2:end)==1;
                    c=[(b(1:end-1)+b(2:end))/2];
                    inshadow2i=[c(1),c,c(end)];
                    kd2(j,:)=kd2(j,:).*(1-inshadow2i);   
                end
                
                %% Plotting
                if 0
                    % plot if needed the situation with the current coastline (in yellow), the structures (black) 
                    % and the crossings that are found with structures (k*) and the coast (ks). 
                    figure(99);clf;plot(COAST.x_mc,COAST.y_mc,'y');hold on;
                    plot(COAST.x,COAST.y,'g');
                    plot(xq(irange),yq(irange),'m.-');   
                    plot(xhard,yhard,'k-');
                    plot([xq(i)+epsxc/1000 xtip(1,j)-epsxc],[yq(i)+epsyc/1000 ytip(1,j)-epsyc],'c-');
                    plot([xq(i)+epsxc/1000 xtip(2,j)-epsxc],[yq(i)+epsyc/1000 ytip(2,j)-epsyc],'c-')
                    plot([xq(i)+epsxc xtip(1,j)-epsxc/1000],[yq(i)+epsyc ytip(1,j)-epsyc/1000],'r-');
                    plot([xq(i)+epsxc xtip(2,j)-epsxc/1000],[yq(i)+epsyc ytip(2,j)-epsyc/1000],'r-');
                    try;plot(xxhard1,yyhard1,'k*');end
                    try;plot(xxhard2,yyhard2,'k*');end
                    try;plot(xx_coast1,yy_coast1,'ks');end
                    try;plot(xx_coast2,yy_coast2,'ks');end
                end
            end
            
            %% reduce impact of very oblique incoming conditions for offshore breakwaters (not for groynes)
            omegamin=45;
            omegamax=80;
            kd1scale = max(1+min((omegamin+om1)/(omegamax-omegamin),0),0);
            kd1scaleR = max(1+min((omegamin-om1)/(omegamax-omegamin),0),0);
            kd2scale = max(1+min((omegamin-om2)/(omegamax-omegamin),0),0);
            kd2scaleR = max(1+min((omegamin+om2)/(omegamax-omegamin),0),0);
            kd1obw=kd1.*kd1scale.*kd2scale.*kd1scaleR;
            kd2obw=kd2.*kd1scale.*kd2scale.*kd2scaleR;
    
            % kdscale=repmat(max(sum((kd1obw.^1+kd2obw.^2).^0.5,1),1),[size(kd1,1),1]);
            % kd1obw=kd1obw./kdscale;
            % kd2obw=kd2obw./kdscale;
    
            %% sumin is the 1d array denoting the influence areas of structures:
            % 0 = no structure, 1 = in lee of one structure, 2 = between two structures
            % and influenced by both.
            in2=in;
            for j=1:nhard
                if STRUC.idgroyne(j)==0
                    in2(j,kd1scaleR(j,:).*kd2scaleR(j,:)==0)=nan;
                end
            end
            sumin=sum(~isnan(in2),1);
            
            %% Points in the overlapping influence areas of two structures
            % multiply kd 
            % do weighted averaging of directions based on wave energy
            om=zeros(size(HStdp));
            irange=find(sumin>=1);
            j_hard_diff=[];
            if ~isempty(irange)
                pwr=2; % scaling of the components of the waves from both sides (by default it is 2 because of wave energy)
                pwr2=2; % used for combining wave heights of different structures (default was always 1, but 2 is more logical value)
                for i=irange
                    j_hard_diff=find(~isnan(in(:,i)));
                    kdj=zeros(size(j_hard_diff));
                    sinkd2=0;
                    coskd2=0;
                    
                    for kk=1:length(j_hard_diff)
                        j=j_hard_diff(kk);
                        if STRUC.idgroyne(j)==0
                            % treat the revetment as an offshore breakwater, but set the component that goes behind the structure to 0.
                            if revetment(j)==1
                                om1(j,itip(1,j)+1:end)=0;
                                om2(j,1:itip(2,j)-1)=0;
                                kd1obw(j,itip(1,j)+1:end)=0;
                                kd2obw(j,1:itip(2,j)-1)=0;
                            end
    
                            % offshore breakwater : add all components
                            kdj(kk)=(kd1obw(j,i).^pwr2+kd2obw(j,i).^pwr2+kd3(j,i).^pwr2);
                            sinkd2=sinkd2+(kd1obw(j,i)*Htip(1,j)).^pwr.*sind(om1(j,i)) ...
                                         +(kd2obw(j,i)*Htip(2,j)).^pwr.*sind(om2(j,i)) ...
                                         +(kd3(j,i)*Htip(3,j)).^pwr.*sind(om3(j,i));
                            coskd2=coskd2+(kd1obw(j,i)*Htip(1,j)).^pwr.*cosd(om1(j,i)) ...
                                         +(kd2obw(j,i)*Htip(2,j)).^pwr.*cosd(om2(j,i)) ...
                                         +(kd3(j,i)*Htip(3,j)).^pwr.*cosd(om3(j,i));
                        elseif (kd1(j,i)*Htip(1,j))~=0
                            % groyne : add left component
                            kdj(kk)=kd1(j,i).^pwr2;
                            sinkd2=sinkd2+(kd1(j,i)*Htip(1,j)).^pwr.*sind(om1(j,i));
                            coskd2=coskd2+(kd1(j,i)*Htip(1,j)).^pwr.*cosd(om1(j,i));
                        else %if (kd2(j,i)*Htip(2,j))~=0
                            % groyne : add right component
                            kdj(kk)=kd2(j,i).^pwr2;
                            sinkd2=sinkd2+(kd2(j,i)*Htip(2,j)).^pwr.*sind(om2(j,i));
                            coskd2=coskd2+(kd2(j,i)*Htip(2,j)).^pwr.*cosd(om2(j,i));
                        end    
                    end 
                    nrhard=1; %max(length(j_hard_diff(j_hard_diff~=0)),1);
                    HStdp(i)=((sinkd2+coskd2)).^(1/pwr2);
                    kdtdp(i)=(sum(kdj)/nrhard).^(1/pwr2);
                    kdtdp2(i)=HStdp(i)/HStdp0(i); % this is a slightly better estimate of kd behind the structure, but has the drawback that it depends on Htip, and therefore have a small step at the end of the diffraction influence area.
                    om(i)=atan2d(sinkd2,coskd2);
                end
            end
            
            %% Reduce gradients at last point of diffraction area
            ids=[find(diff(sumin>=1)==1)+1,find(diff(sumin>=1)==-1)];
            ids=setdiff(ids,[2,length(kdtdp)-1]);
            kdtdp(ids)=(1+kdtdp(ids))/2;
            om(ids)=om(ids)/2;
            
            %% Smooth diffraction area
            om=WAVE.diffsmooth*get_smoothdata(om,'')+(1-WAVE.diffsmooth)*om;
            om(~sumin)=0;
            kdtdp(~sumin)=1;
    
            %% OUTPUT:
            %WAVE.dirtip(qq,:)   = dirtip(:)';
            WAVE.diff(qq,:)   = sumin;                        % wave diffraction at QS points
            % WAVE.diffxy(qq,:) = [sumin(:,1),sumin(:,1:end-1)&sumin(:,2:end),sumin(:,end)];  % wave diffraction at XY points
            % WAVE.diff2{qq}    = in;
            % STRUC.xtip(qq,:)  = xtip(:)'; 
            % STRUC.ytip(qq,:)  = ytip(:)'; 
            % STRUC.Htip(qq,:)  = Htip(:)'; 
            STRUC.xp          = xq; 
            STRUC.yp          = yq; 
            WAVE.PHItdp(qq,:) = PHItdp-om;  % minus because om is cartesian, PHItdp nautical
            WAVE.HStdp(qq,:)  = HStdp0.*min(kdtdp,1);
            
            %% Overview plot of impact of diffraction on wave conditions (relative wave height modification and wave angle change)
            % this is the diagnostic plot for wave diffraction
            if 0
                figure(24+COAST.i_mc);clf;
                ha1=subplot(2,1,1);
                nhardj=length(j_hard_diff);
                col=jet(nhard);
                lnst={'-','-','-','-'};
                sq=[0,cumsum([(diff(COAST.xq).^2+diff(COAST.yq).^2).^0.5])];
                hp=plot(sq/1000,kdtdp,'g-','LineWidth',1);hold on;
                legtxt={'combined'};
                for jj0=1:nhard
                    jj=jj0;
                    hp(jj0+1)=plot(sq/1000,kd1plot(jj,:),'Color',col(jj0,:)*0.5,'LineStyle',lnst{mod(jj0-1,4)+1});
                    plot(sq/1000,kd2plot(jj,:),'Color',col(jj0,:)*0.5,'LineStyle',lnst{mod(jj0-1,4)+1});
                    legtxt{jj0+1}=['Structure ',num2str(jj)];
                end
                plot(sq/1000,kdtdp,'g:','LineWidth',2);hold on;
                %plot(sq/1000,kdtdp2,'g:','LineWidth',1.5);hold on;
                ylabel('kd value [-]');
                xlim([min(sq),max(sq)]/1000);
                pos1=get(gca,'Position');
                ha2=subplot(2,1,2);
                hp=plot(sq/1000,om,'g-','LineWidth',1);hold on;
                for jj0=1:nhard
                    jj=jj0;
                    hp(jj0+1)=plot(sq/1000,om1plot(jj,:),'Color',col(jj0,:)*0.5,'LineStyle',lnst{mod(jj0-1,4)+1});hold on;
                    plot(sq/1000,om2plot(jj,:),'Color',col(jj0,:)*0.5,'LineStyle',lnst{mod(jj0-1,4)+1});
                end
                plot(sq/1000,om,'g:','LineWidth',2);hold on;
                xlabel('alongshore distance [km]');
                ylabel('wave dir rotation [°]');
                xlim([min(sq),max(sq)]/1000);
                pos2=get(gca,'Position');
                legend(hp,legtxt,'Location','SouthOutside','Orientation','Horizontal');
                warning off
                set(ha1,'Position',[0.08,0.58,0.85,0.37]);
                set(ha2,'Position',[0.08,0.16,0.85,0.37]);
                set(gcf,'Position',[10+(COAST.i_mc-1)*1900/COAST.n_mc,600,1900/COAST.n_mc,400]);
                warning on
                figure(11)
            end      
        end
    end
end

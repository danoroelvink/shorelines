function [V,FORMAT,TIME]=plot_coast(CHANNEL,STRUC,COAST,DUNE,WAVE,TIME,TRANSP,FORMAT,V,FNOUR)
% function [V,FORMAT,TIME]=plot_coast(CHANNEL,STRUC,COAST,DUNE,WAVE,TIME,TRANSP,FORMAT,V,FNOUR)
%
% INPUT:
%   STRUC
%      .x_hard
%      .y_hard
%   COAST
%      .x_mc0
%      .y_mc0
%      .PHIf_mc
%      .x_mc
%      .y_mc
%      .Wberm_mc
%   WAVE
%      .HSo_mc
%      .PHIo_mc
%      .HSbr_mc
%      .PHIbr_mc
%      .WVC
%   TIME
%      .it
%      .tnow
%      .dt
%      .tnext
%   TRANSP
%      .QS_mc
%   FORMAT
%      .tplot
%      .tsl
%      .plot_time
%      .xmax
%      .ymax
%      .xmin
%      .ymin
%      .nmax
%      .xp_mc
%      .yp_mc
%      .xlimits
%      .ylimits
%      .XYoffset
%      .SLplot
%      .fignryear
%      .outputdir
%      .video
%      .plotinterval;
%      .usefill;
%      .LDBplot;
%      .llocation
%      .ld
%   V
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
%   This library is free software: you can redistribute TIME.it and/or
%   modify TIME.it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation, either
%   version 2.1 of the License, or (at your option) any later version.
%
%   This library is distributed in the hope that TIME.it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses
%   --------------------------------------------------------------------

    qvrscale=1.5;
    eps=-1d-10;
    cols=[0.1 0.1 0.1; 0 0.6 0;0.5 0.1 0;0 0.2 0.5]; % colours of respectively the waves offshore (OFF), at the depth-of-closure (TDP) and at the point of breaking (BR)    
    if isempty(WAVE.HStdp_mc)
        WAVE.HStdp_mc=WAVE.HSo_mc.*0;
        WAVE.PHIbr_mc=WAVE.PHIo_mc;
    end
    if isempty(WAVE.HSbr_mc)
        WAVE.HSbr_mc=WAVE.HSo_mc.*0;
        WAVE.PHIbr_mc=WAVE.PHIo_mc;
    end
    WAVE.HSbr_mc(isnan(WAVE.HSbr_mc))=0;
    IDnotnan=~isnan(WAVE.HSo_mc);
    HS_mean = median(WAVE.HSo_mc(WAVE.HSo_mc-mean(WAVE.HSo_mc(IDnotnan))>eps));       % necessary for spatially uniform waves
    PHI_mean = mod(atan2d(median(sind(WAVE.PHIo_mc(IDnotnan))),median(cosd(WAVE.PHIo_mc(IDnotnan)))),360);
    HStdp_mean = median(WAVE.HStdp_mc(WAVE.HSo_mc-mean(WAVE.HSo_mc(IDnotnan))>eps));
    PHItdp_mean = mod(atan2d(median(sind(WAVE.PHItdp_mc(IDnotnan))),median(cosd(WAVE.PHItdp_mc(IDnotnan)))),360);
    HSbr_mean = median(WAVE.HSbr_mc(WAVE.HSo_mc-mean(WAVE.HSo_mc(IDnotnan))>eps));
    PHIbr_mean = mod(atan2d(median(sind(WAVE.PHIbr_mc(IDnotnan))),median(cosd(WAVE.PHIbr_mc(IDnotnan)))),360);
    shadowS_mc=TRANSP.shadowS_mc;
    shadowS_mc(isnan(shadowS_mc))=1;
    shadowS_h_mc=TRANSP.shadowS_h_mc;
    shadowS_h_mc(isnan(shadowS_h_mc))=1;
    if length(shadowS_mc)~=length(WAVE.HStdp_mc)
        shadowS_mc=zeros(1,length(WAVE.HStdp_mc));
        shadowS_h_mc=zeros(1,length(WAVE.HStdp_mc));
    end
    if ~exist(fullfile(pwd,FORMAT.outputdir),'dir')
        mkdir(fullfile(pwd,FORMAT.outputdir));
    end

    if mod(TIME.it,FORMAT.plotinterval)==0
        %iplot=iplot+1;
        if 1
            %% BACKGROUND LAYER
            try
                hold off
                [h,ax] = showgeorefimage(FORMAT.IMGfilename,FORMAT.WORLDfilename);hold on;
                xlim(FORMAT.xlimits);
                ylim(FORMAT.ylimits);
            end
            
            %% FILL THE COASTLINE
            [FORMAT] = plot_fill_sections(COAST,FORMAT,TIME);
            hold on;
            
            %% Plot shoreface nourishment 
            if FNOUR.fnourish == 1	            
                map = colormap(jet(FNOUR.n));
                for i = 1 : FNOUR.n
                    plot( FNOUR.x(i,:) , FNOUR.y(i,:) , 'color' , map(i,:) , 'linewidth' , 2 ); 	                    
                end 
            end 
            
            %% PLOT COASTAL STRUCTURES
            plot(COAST.x_mc0,COAST.y_mc0,'k-','linewidth',1,'Color',[0.5 0.5 0.5]); 
            plot(COAST.x_mc,COAST.y_mc,'b-','linewidth',1,'Color',[0.2 0.2 0.8]); 
            for ist=1:STRUC.n_hard
                [xhard,yhard] = get_one_polygon( STRUC.x_hard,STRUC.y_hard,ist );
                if STRUC.idgroyne(ist)~=0
                    hf=fill(xhard,yhard,[0.5 0.5 0.3]);
                    set(hf,'LineStyle','none');
                elseif (abs(xhard(end)-xhard(1))+abs(yhard(end)-yhard(1)))<1e-2
                    hf=fill(xhard,yhard,[0.5 0.5 0.3]);
                end
                plot(xhard,yhard,'k-','linewidth',1);%try id1=[1,find(isnan(STRUC.x_hard))+1]; plot(STRUC.x_hard(id1),STRUC.y_hard(id1),'ks');end
            end
            plot(STRUC.x_revet,STRUC.y_revet,'m-','linewidth',2);%try id1=[1,find(isnan(STRUC.x_revet))+1]; plot(STRUC.x_revet(id1),STRUC.y_revet(id1),'ms');end
            if TRANSP.sedlim
                plot(TRANSP.x_sedlim,TRANSP.y_sedlim,'c:','linewidth',2);%try id1=[1,find(isnan(STRUC.x_revet))+1]; plot(STRUC.x_revet(id1),STRUC.y_revet(id1),'ms');end
            end
            if ~isempty(CHANNEL.xr_mc') 
                plot(CHANNEL.xr_mc,CHANNEL.yr_mc,'b.-')
                plot(CHANNEL.x_inlet,CHANNEL.y_inlet,'*')
            end

            %% PLOT DUNES
            if DUNE.used==1
                xdune0 = DUNE.xdune0; 
                ydune0 = DUNE.ydune0; 
                xdune = COAST.x_mc - COAST.Wberm_mc.*sind(COAST.PHIcxy_mc); 
                ydune = COAST.y_mc - COAST.Wberm_mc.*cosd(COAST.PHIcxy_mc); 
                hpd0=plot(xdune0,ydune0,'k--','linewidth',1,'Color',[0.5 0.5 0.5]);
                hpd=plot(xdune,ydune,'k--','linewidth',1,'Color',[0.2 0.2 0.8]);
            end

            %% PLOT OFFSHORE WAVE HEIGHT TEXT
            if isempty(FORMAT.XYwave)
                FORMAT.XYwave = [FORMAT.XYoffset(1)+FORMAT.xlimits(1)+(sind(PHI_mean)+1)/2*diff(FORMAT.xlimits),FORMAT.XYoffset(2)+FORMAT.ylimits(1)+(cosd(PHI_mean)+1)/2*0.95*diff(FORMAT.ylimits),diff(FORMAT.xlimits)/12];
            end

            %% PLOT WAVE QUIVER & OFFSHORE WAVE HEIGHT TEXT
            PHI2=[];
            PHItdp2=[];
            PHIbr2=[];
            PHIf2=[];
            HS2=[];
            HStdp2=[];
            HSbr2=[];
            xx=[];
            yy=[];  
            if isempty(WAVE.WVC) || WAVE.spacevaryingwave==0
                arx=[FORMAT.XYwave(1)-FORMAT.XYoffset(1)];
                ary=[FORMAT.XYwave(2)-FORMAT.XYoffset(2)];
                PHI2=PHI_mean;
                PHItdp2=PHItdp_mean;
                PHIbr2=PHIbr_mean;
                HS2=HS_mean;
                HStdp2=HStdp_mean;
                HSbr2=HSbr_mean;
                if ~isempty(COAST.PHIf_mc) 
                    PHIf2=median(COAST.PHIf_mc(~isnan(COAST.PHIf_mc)));
                end
                xx=min(max(FORMAT.XYwave(1)-FORMAT.XYoffset(1)-diff(FORMAT.xlimits)/24,min(FORMAT.xlimits)),max(FORMAT.xlimits));
                yy=min(max(FORMAT.XYwave(2)-FORMAT.XYoffset(2)+diff(FORMAT.ylimits)/24,min(FORMAT.ylimits)),max(FORMAT.ylimits));
            else %if ~isempty(COAST.x_mc)
                WVCid=zeros(length(WAVE.WVC),1);
                arx=zeros(length(WAVE.WVC),1);
                ary=zeros(length(WAVE.WVC),1);
                PHI2=zeros(length(WAVE.WVC),1);
                PHItdp2=zeros(length(WAVE.WVC),1);
                PHIbr2=zeros(length(WAVE.WVC),1);
                HS2=zeros(length(WAVE.WVC),1);
                HStdp2=zeros(length(WAVE.WVC),1);
                HSbr2=zeros(length(WAVE.WVC),1);
                xx=zeros(length(WAVE.WVC),1);
                yy=zeros(length(WAVE.WVC),1);
                PHIf2=zeros(length(WAVE.WVC),1);
                for gg=1:length(WAVE.WVC)
                    try
                        dist=((COAST.xq1_mc-WAVE.WVC(gg).x).^2+(COAST.yq1_mc-WAVE.WVC(gg).y).^2).^0.5;
                        WVCid(gg)=min(find(dist==min(dist),1),length(WAVE.PHIo_mc));
                        arx(gg)=[WAVE.WVC(gg).x-FORMAT.XYoffset(1)];
                        ary(gg)=[WAVE.WVC(gg).y-FORMAT.XYoffset(2)];
                        PHI2(gg)=WAVE.PHIo_mc(WVCid(gg));
                        PHItdp2(gg)=WAVE.PHItdp_mc(WVCid(gg));
                        PHIbr2(gg)=WAVE.PHIbr_mc(WVCid(gg));
                        HS2(gg)=WAVE.HSo_mc(WVCid(gg));
                        HStdp2(gg)=WAVE.HStdp_mc(WVCid(gg));
                        HSbr2(gg)=WAVE.HSbr_mc(WVCid(gg));
                        xx(gg)=arx(gg);
                        yy(gg)=ary(gg);
                        if ~isempty(COAST.PHIf_mc)
                            PHIf2(gg)=COAST.PHIf_mc(min(WVCid(gg),length(COAST.PHIf_mc)));
                            xx(gg)=arx(gg) + diff(FORMAT.xlimits)/35 .* sind(PHIf2(gg));
                            yy(gg)=ary(gg) + diff(FORMAT.ylimits)/35 .* cosd(PHIf2(gg));
                        end
                    end
                end
            end
            for gg=1:length(PHI2)
                if arx(gg)>FORMAT.xlimits(1) && arx(gg)<FORMAT.xlimits(2) && ary(gg)>FORMAT.ylimits(1) && ary(gg)<FORMAT.ylimits(2) 
                    factor = HS2(gg)/max(HS2);
                    hqvr1=quiver(arx(gg),ary(gg),FORMAT.XYwave(3).*cosd(3/2*180-PHI2(gg))*factor,FORMAT.XYwave(3).*sind(3/2*180-PHI2(gg))*factor,qvrscale);hold on;
                    set(hqvr1,'linewidth',2,'Color',cols(1,:),'AutoScale','off','AutoScaleFactor',1.2,'MaxHeadSize',1);
                    hqvr2=quiver(arx(gg),ary(gg),FORMAT.XYwave(3).*cosd(3/2*180-PHItdp2(gg))*HStdp2(gg)/HS2(gg)*factor,FORMAT.XYwave(3).*sind(3/2*180-PHItdp2(gg))*HStdp2(gg)/HS2(gg)*factor,qvrscale);hold on;
                    set(hqvr2,'linewidth',2,'Color',cols(2,:),'AutoScale','off','AutoScaleFactor',1.2,'MaxHeadSize',1);
                    hqvr3=quiver(arx(gg),ary(gg),FORMAT.XYwave(3).*cosd(3/2*180-PHIbr2(gg))*HSbr2(gg)/HS2(gg)*factor,FORMAT.XYwave(3).*sind(3/2*180-PHIbr2(gg))*HSbr2(gg)/HS2(gg)*factor,qvrscale);hold on;
                    set(hqvr3,'linewidth',2,'Color',cols(3,:),'AutoScale','off','AutoScaleFactor',1.2,'MaxHeadSize',1);
                    htxt=text(xx(gg)+FORMAT.XYwave(3).*cosd(PHIf2(gg))/3,yy(gg)+FORMAT.XYwave(3).*sind(PHIf2(gg))/3,[' ',num2str(HSbr2(gg),'%2.1f'),'m']);
                    set(htxt,'FontSize',9,'HorizontalAlignment','Center','VerticalAlignment','Middle','Color',cols(1,:));
                end
            end

            %% PLOT LOWER SHOREFACE ORIENTATION
            if ~isempty(PHIf2) && ~strcmpi(TRANSP.trform,'KAMP') && ~strcmpi(TRANSP.trform,'CERC') && ~strcmpi(TRANSP.trform,'CERC2')
                for gg=1:length(arx)
                    dx=FORMAT.XYwave(3).*cosd(PHIf2(gg))/4.*[1,-1];
                    dy=FORMAT.XYwave(3).*sind(PHIf2(gg))/4.*[-1,1];
                    hl=plot(arx(gg)+dx,ary(gg)+dy);hold on;
                    set(hl,'linewidth',1,'Color',[0 0.5 0]);
                end
            end
            
            %% PLOT WAVE HEIGHT AT DEPTH-OF-CLOSURE (optional)
            if FORMAT.plotHS>0 && ~isempty(WAVE.HStdp_mc)
                hold on;
                hsc1=scatter(COAST.xq1_mc-FORMAT.XYoffset(1),COAST.yq1_mc-FORMAT.XYoffset(2),6,WAVE.HSbr_mc);
                set(hsc1,'markerFaceColor','flat');
                set(gca,'Clim',[0,5]);hold on;
                ggrange=[1:FORMAT.plotHS:min(length(WAVE.HStdp_mc),length(COAST.xq1_mc))];
                hhs=plot(COAST.xq1_mc(ggrange),COAST.yq1_mc(ggrange),'kx');
                for gg=ggrange
                    idxy=COAST.xq1_mc(gg)>FORMAT.xlimits(1) && COAST.xq1_mc(gg)<FORMAT.xlimits(2) && COAST.yq1_mc(gg)>FORMAT.ylimits(1) && COAST.yq1_mc(gg)<FORMAT.ylimits(2);
                    if idxy && ~isnan(COAST.xq1_mc(gg)) && shadowS_h_mc(gg)<1 && shadowS_mc(gg)<1
                        htxt=text(COAST.xq1_mc(gg),COAST.yq1_mc(gg)-FORMAT.XYoffset(2),['',num2str(WAVE.HStdp_mc(gg),'%2.1f'),'m']);
                        set(htxt,'FontSize',9,'HorizontalAlignment','Center','VerticalAlignment','Bottom','Color',cols(2,:));
                        htxt=text(COAST.xq1_mc(gg),COAST.yq1_mc(gg)-FORMAT.XYoffset(2),['',num2str(WAVE.HSbr_mc(gg),'%2.1f'),'m']);
                        set(htxt,'FontSize',9,'HorizontalAlignment','Center','VerticalAlignment','Top','Color',cols(3,:));
                    end
                end
            end
            
            %% PLOT WAVE DIRECTION AT DEPTH-OF-CLOSURE (optional)
            if FORMAT.plotDIR>0 && ~isempty(WAVE.PHItdp_mc)
                hold on;                
                for gg=[max(round(FORMAT.plotDIR*0.6),1):FORMAT.plotDIR:min(length(TRANSP.QS_mc),length(COAST.xq1_mc))]
                    for mm=1:2                       
                        vertalignm={'Bottom','Top'};
                        if mm==1
                            waveheight=WAVE.HStdp_mc(gg);
                            PHIwave=WAVE.PHItdp_mc(gg);
                        else
                            waveheight=WAVE.HSbr_mc(gg);
                            PHIwave=WAVE.PHIbr_mc(gg);
                        end   
                        idxy=COAST.xq1_mc(gg)>FORMAT.xlimits(1) && COAST.xq1_mc(gg)<FORMAT.xlimits(2) && COAST.yq1_mc(gg)>FORMAT.ylimits(1) && COAST.yq1_mc(gg)<FORMAT.ylimits(2);
                        if idxy && ~isnan(COAST.xq1_mc(gg)) && shadowS_h_mc(gg)<1 && shadowS_mc(gg)<1
                            factor = sqrt(waveheight/2)*0.7;
                            QVR1=diff(FORMAT.xlimits)/12.*cosd(3/2*180-PHIwave).*factor;
                            QVR2=diff(FORMAT.xlimits)/12.*sind(3/2*180-PHIwave).*factor;
                            hqvr9=quiver(COAST.xq1_mc(gg)-FORMAT.XYoffset(1),COAST.yq1_mc(gg)-FORMAT.XYoffset(2),QVR1,QVR2,(max(1-shadowS_h_mc(gg)-shadowS_mc(gg),0))*qvrscale);hold on;
                            set(hqvr9,'linewidth',2,'Color',cols(mm+1,:),'AutoScale','off','AutoScaleFactor',1.2,'MaxHeadSize',1);               
                            htxt=text(COAST.xq1_mc(gg)-FORMAT.XYoffset(1),COAST.yq1_mc(gg)-FORMAT.XYoffset(2),['',num2str(PHIwave,'%1.0f'),'°']);
                            set(htxt,'FontSize',9,'HorizontalAlignment','Center','VerticalAlignment',vertalignm{mm},'Color',cols(mm+1,:));
                        end
                    end
                end
            end
            
            %% PLOT TRANSPORTS (optional)
            if FORMAT.plotQS>0 && ~isempty(TRANSP.QS_mc)
                hold on;
                hsc1=scatter(COAST.xq1_mc-FORMAT.XYoffset(1),COAST.yq1_mc-FORMAT.XYoffset(2),6,TRANSP.QS_mc);
                set(hsc1,'markerFaceColor','flat');
                set(gca,'Clim',[-1,1]*200000);hold on;
                ggrange=[1:FORMAT.plotQS:min(length(TRANSP.QS_mc),length(COAST.xq1_mc))];
                hqs=plot(COAST.xq1_mc(ggrange),COAST.yq1_mc(ggrange),'k+');
                for gg=ggrange
                    idxy=COAST.xq1_mc(gg)>FORMAT.xlimits(1) && COAST.xq1_mc(gg)<FORMAT.xlimits(2) && COAST.yq1_mc(gg)>FORMAT.ylimits(1) && COAST.yq1_mc(gg)<FORMAT.ylimits(2);
                    if idxy && ~isnan(COAST.xq1_mc(gg))
                        htxt=text(COAST.xq1_mc(gg)-FORMAT.XYoffset(1),COAST.yq1_mc(gg)-FORMAT.XYoffset(2),['',num2str(TRANSP.QS_mc(gg)/1000,'%1.0f'),'m3/yr']);
                        set(htxt,'FontSize',9,'HorizontalAlignment','Center','VerticalAlignment','Middle','Color',cols(1,:));
                    end
                end
            end
            
            %% PLOT DIFFRACTION POINTS
            hpd=plot(COAST.xq1_mc(WAVE.diff_mc)-FORMAT.XYoffset(1),COAST.yq1_mc(WAVE.diff_mc)-FORMAT.XYoffset(2),'.','MarkerSize',8,'Color',[0.5 0.5 0.3]);

            %% PLOT UPWIND CORRECTION POINTS
            if FORMAT.plotUPW>0 && ~isempty(TRANSP.ivals_mc)
                if length(TRANSP.ivals_mc)==length(COAST.xq1_mc)
                    hup1=plot(COAST.xq1_mc(TRANSP.ivals_mc==1),COAST.yq1_mc(TRANSP.ivals_mc==1),'rs');
                    hup2=plot(COAST.xq1_mc(TRANSP.ivals_mc==-1),COAST.yq1_mc(TRANSP.ivals_mc==-1),'gs');
                end
            end
            
            %% PLOT REFERENCE LINES AND LEGEND
            if ~isempty(FORMAT.LDBplot)
                for mm=1:size(FORMAT.LDBplot,1)
                    try
                        LDBplotval = get_landboundary(FORMAT.LDBplot{mm,1});
                    catch
                        LDBplotval = load(FORMAT.LDBplot{mm,1});
                    end
                    hp6(mm)=plot(LDBplotval(:,1)-FORMAT.XYoffset(1),LDBplotval(:,2)-FORMAT.XYoffset(2),FORMAT.LDBplot{mm,3},'linewidth',0.5);
                end
                hleg = legend(hp6,FORMAT.LDBplot(:,2)','Location',FORMAT.llocation);
                set(hleg,'Box','off','Color','None');
            end
            
            %% PLOT PROFILES
            if ~isempty(FORMAT.xyprofiles)
                plot(FORMAT.xyprofiles(:,1),FORMAT.xyprofiles(:,2),'gs');
                plot(FORMAT.xyprofiles(:,1),FORMAT.xyprofiles(:,2),'g+');
                for pp=1:size(FORMAT.xyprofiles,1)
                    dist=((COAST.x_mc-FORMAT.xyprofiles(pp,1)).^2+(COAST.y_mc-FORMAT.xyprofiles(pp,2)).^2).^0.5;
                    idp=find(dist==min(dist),1);
                    xduneprofile = [COAST.x_mc(idp), COAST.x_mc(idp) - COAST.Wberm_mc(idp).*sind(COAST.PHIcxy_mc(idp))]; 
                    yduneprofile = [COAST.y_mc(idp), COAST.y_mc(idp) - COAST.Wberm_mc(idp).*cosd(COAST.PHIcxy_mc(idp))]; 
                    plot(xduneprofile,yduneprofile,'g--');
                    htxt=text(FORMAT.xyprofiles(pp,1),FORMAT.xyprofiles(pp,2),['profile',num2str(pp)]);
                    set(htxt,'FontSize',9,'HorizontalAlignment','Right','VerticalAlignment','Top');
                end
            end
            
            %% PLOT SHORELINES AT SPECIFIC DATES
            if ~isempty(FORMAT.SLplot)&& TIME.dt==(FORMAT.tplot-TIME.tnow)/365
                try
                FORMAT.xp_mc{FORMAT.tsl,:}=COAST.x_mc;
                FORMAT.yp_mc{FORMAT.tsl,:}=COAST.y_mc;
                FORMAT.tsl=FORMAT.tsl+1;
                FORMAT.tplot=FORMAT.plot_time(FORMAT.tsl);
                end
            else%if isempty(FORMAT.SLplot)
                FORMAT.xp_mc=FORMAT.xp_mc;
                FORMAT.yp_mc=FORMAT.yp_mc;
            end
            if ~isempty(FORMAT.SLplot) && FORMAT.tsl>1
                if ~isempty(FORMAT.LDBplot)
                    LName=FORMAT.LDBplot(:,2);
                else
                    LName={};
                    mm=0;
                end
                for sl=1:size(FORMAT.xp_mc,1)
                    hp6(sl+mm)=plot(FORMAT.xp_mc{sl,:},FORMAT.yp_mc{sl,:},FORMAT.SLplot{sl,3},'linewidth',1.5);
                    LName(end+1)=FORMAT.SLplot(sl,2);
                end
                hleg = legend(hp6,LName','Location',FORMAT.llocation);
                set(hleg,'Box','off','Color','None');
            elseif isempty(FORMAT.SLplot) && ~isempty(FORMAT.LDBplot)
                hleg = legend(hp6,FORMAT.LDBplot(:,2)','Location',FORMAT.llocation);
                set(hleg,'Box','off','Color','None');
            end
        end
        
        %% FORMATTING THE PLOT
        hold off;
        ax = get(FORMAT.mainfighandle,'CurrentAxes');
        set(ax,'DataAspectRatio',[1 1 1],'DataAspectRatioMode','manual','PlotBoxAspectRatio',[3 4 4],'PlotBoxAspectRatioMode','manual'); % manual axis equal
        xlim(FORMAT.xlimits);
        ylim(FORMAT.ylimits);
        if max(xlim)>=600
        set(ax,'XtickLabel',num2str(get(ax,'Xtick')'/1000,'%2.1f'));
        else
        set(ax,'XtickLabel',num2str(get(ax,'Xtick')'/1000,'%2.2f'));
        end    
        if max(ylim)>=600
        set(ax,'YtickLabel',num2str(get(ax,'Ytick')'/1000,'%2.1f'));
        else
        set(ax,'YtickLabel',num2str(get(ax,'Ytick')'/1000,'%2.2f'));
        end
        xlabel(ax,'Easting [km]','HandleVisibility','off');
        ylabel(ax,'Northing [km]','HandleVisibility','off');
        date=datevec(double(TIME.tnow));
        title(ax,num2str(date(1:3)));
        if ~FORMAT.fastplot
           drawnow;    % 10% of plottime
        end
        
        %% video
        if FORMAT.video==1 && STRUC.diffraction==0
            V(TIME.it+1)=getframe(FORMAT.mainfighandle);
        end
        
        %% images
        if TIME.tnow>=TIME.tnext
            fname=[num2str(round((TIME.it+1)),'%05.0f')]; 
            if ~FORMAT.fastplot    
               print(FORMAT.mainfighandle,fullfile(FORMAT.outputdir,[fname '.png']),'-dpng','-r150');
            else
               F=getframe(gcf);
               imwrite(F.cdata,fullfile(FORMAT.outputdir,[fname '.png']),'png');
            end
            vii=length(V)+1;
            V(vii)=getframe(FORMAT.mainfighandle);
            TIME.tnext=TIME.tnext+365./FORMAT.fignryear;
        end
        FORMAT.xyt(TIME.it+1).COAST.x_mc=COAST.x_mc;
        FORMAT.xyt(TIME.it+1).COAST.y_mc=COAST.y_mc;
        
        drawnow;
    end
end

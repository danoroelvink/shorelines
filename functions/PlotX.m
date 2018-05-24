%function [vi,xp_mc,yp_mc,S,nmax,tsl,tplot]=PlotX(S,x_hard,y_hard,x_mc0,y_mc0,it,x_mc,y_mc,ld,tnow,phiw,tplot,tsl,plot_time,nmax,xp_mc,yp_mc,vi)
% function [vi,xp_mc,yp_mc,S,xmax,ymax,xmin,ymin,nmax,tsl,tplot,vii]=PlotX(S,x_hard,y_hard,x_mc0,y_mc0,it,x_mc,y_mc,ld,tnow,phiw,tplot,tsl,plot_time,xmax,ymax,xmin,ymin,nmax,xp_mc,yp_mc,vi,vii,xp,yp,QS,xS)
function [vi,xp_mc,yp_mc,S,xmax,ymax,xmin,ymin,nmax,tsl,tplot,vii]=PlotX(S,x_hard,y_hard,x_mc0,y_mc0,it,x_mc,y_mc,ld,tnow,phiw,tplot,tsl,plot_time,xmax,ymax,xmin,ymin,nmax,xp_mc,yp_mc,vi,vii)

if mod(it,S.plotinterval)==0
    %iplot=iplot+1;
    if 1
%         subplot(211)
         [xmax,ymax,xmin,ymin,nmax ] = fill_sections(x_mc,y_mc,ld,it,xmax,ymax,xmin,ymin,nmax);
       % [nmax] = fill_sections(x_mc,y_mc,ld,it,nmax);
        
        %axis([-10000 10000 0 10000])
        %set(gca,'xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[])
        hold on;
        plot(x_hard,y_hard,'k','linewidth',2);
        plot(x_mc0,y_mc0,'b','linewidth',2); %2
%         plot(xp,yp,'*');
        if isempty(S.XYwave)
            S.XYwave = [S.xlimits(1)+(sin(S.phiw0)+1)/2*diff(S.xlimits),S.ylimits(1)+(cos(S.phiw0)+1)/2*diff(S.ylimits),diff(S.xlimits)/12];
        end
        if ~isempty(S.reftime)
            text(S.XYwave(1)-S.XYoffset(1),S.XYwave(2)-S.XYoffset(2)+diff(S.ylimits)/12,['H_{m0}=',num2str(max(S.Hso),'%2.1f'),' (',datestr(tnow,'yyyy-mmm'),')']);
        else
            text(S.XYwave(1)-S.XYoffset(1),S.XYwave(2)-S.XYoffset(2)+diff(S.ylimits)/12,['H_{m0}=',num2str(max(S.Hso),'%2.1f'),' (',num2str(tnow*24*60*60,'%3.2f'),' sec)']);
        end
        if 1
            arx=[S.XYwave(1)-S.XYoffset(1),S.XYwave(1)-S.XYoffset(1)+S.XYwave(3)*cos(3/2*pi-phiw)];
            ary=[S.XYwave(2)-S.XYoffset(2),S.XYwave(2)-S.XYoffset(2)+S.XYwave(3)*sin(3/2*pi-phiw)];
            %plot(arx,ary,'linewidth',2);
            qvrscale=1.5;
            hqvr=quiver(arx(1),ary(1),S.XYwave(3)*cos(3/2*pi-phiw),S.XYwave(3)*sin(3/2*pi-phiw),qvrscale);
            set(hqvr,'linewidth',2,'Color','k','AutoScale','on','AutoSCaleFactor',1.5);
            %set(hqvr,'MaxHeadSize',50,'MarkerSize',100);
        end
        
        %% PLOT REFERENCE LINES AND LEGEND
        if ~isempty(S.LDBplot)
            for mm=1:size(S.LDBplot,1)
                try
                    LDBplotval = landboundary('read',S.LDBplot{mm,1});
                catch
                    LDBplotval = load(S.LDBplot{mm,1});
                end
                hp6(mm)=plot(LDBplotval(:,1)-S.XYoffset(1),LDBplotval(:,2)-S.XYoffset(2),S.LDBplot{mm,3},'linewidth',1.5);%,'MarkerSize',10,'MarkerIndices',1:10:length(x_mc0)); %for test1
            end
%             hleg = legend(hp6,S.LDBplot(:,2)','Location',S.legendlocation);
%             set(hleg,'Box','off','Color','None');
        end
        
        %% Plot shorelines at specific dates
        
        if ~isempty(S.SLplot)&& S.dt==(tplot-tnow)/365
            xp_mc{tsl,:}=x_mc;
            yp_mc{tsl,:}=y_mc;
            tsl=tsl+1;
            tplot=plot_time(tsl);
        else%if isempty(S.SLplot)
            xp_mc=xp_mc;
            yp_mc=yp_mc;
        end
        if ~isempty(S.SLplot) && tsl>1
            if ~isempty(S.LDBplot)
                LName=S.LDBplot(:,2);
            else
                LName={};
                mm=0;
            end
            
            for sl=1:size(xp_mc,1)
                hp6(sl+mm)=plot(xp_mc{sl,:},yp_mc{sl,:},S.SLplot{sl,3},'linewidth',1.5);
                LName(end+1)=S.SLplot(sl,2);
            end
%             hleg = legend(hp6,LName','Location',S.legendlocation);
%             set(hleg,'Box','off','Color','None');
        elseif isempty(S.SLplot) && ~isempty(S.LDBplot)
            hleg = legend(hp6,S.LDBplot(:,2)','Location',S.legendlocation);
            set(hleg,'Box','off','Color','None');
        end
        %for shorelines plotting and extracting
        
        
        %% FORMATTING
        hold off;
        axis equal;
        axis 'auto y'
        xlim(S.xlimits);
        ylim(S.ylimits);
        set(gca,'XtickLabel',num2str(get(gca,'Xtick')'/1000,'%2.1f'));
        set(gca,'YtickLabel',num2str(get(gca,'Ytick')'/1000,'%2.1f'));
%         set(gca,'DataAspectRatio',[0.25 0.0625 1]) %for test1
        xlabel('Easting [km]');
        ylabel('Northing [km]');
%         date=datevec(tnow);
%         title(num2str(date(1:3)));
        drawnow;
        
        %video
        
%         vi(it+1)=getframe(figure(11));
%         pause(0.01);
% subplot(212)
% plot(xS,QS)

        figstep=round((1/S.fignryear)/S.dt);
        if mod(it,figstep)==0
            %fname=[num2str(it+1000),'.jpg'];
            fname=[num2str(round((it+1)+100000)),'.jpg'];
            setFIGUREproperties(800,600,32);
            if ~exist(S.outputdir)
                mkdir(S.outputdir);
            end
            
             print('-djpeg','-r300',[S.outputdir,filesep,fname]);
             vii=vii+1;
             vi(vii)=getframe(figure(11));
             pause(0.01);
        
        end
        %         subplot(122)
        %         plot(philoc(xS<1200),yS(xS<1200),'o')
        %         ylim([2000 8000])
        %         title(['section ',num2str(i_mc)])
        %         %xlim([-1e6 1e6])
              

        
        
    end
    S.xyt(it+1).x_mc=x_mc;
    S.xyt(it+1).y_mc=y_mc;
    
end
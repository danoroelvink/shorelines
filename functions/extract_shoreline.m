function []=extract_shoreline(S,x_hard,y_hard,x_mc0,y_mc0,xp_mc,yp_mc)
if S.print_fig
    if strcmpi(S.trform,'KAMP')
 print((strcat(S.trform,'_',num2str(2.33),'_ds0_',num2str(S.ds0),'_phi_',num2str(S.phiw0/pi*180))),'-dpng', '-r300')
%   print((strcat('xhrd_',num2str(S.x_hard),'_yhrd_',num2str(S.y_hard),'_phi0_2',num2str(S.phiw0))),'-dpng', '-r300')
    else  
print((strcat(S.trform,'_',num2str(S.b/1e4),'_ds0_',num2str(S.ds0),'_phi_',num2str(S.phiw0/pi*180))),'-dpng', '-r300')

    end
end

if ~isempty(S.SLplot)
    if ~isempty(S.LDBplot)
        L2Name(:)=S.LDBplot(:,2);
    else
        L2Name(:)={};
        mm=0;
    end
    figure(12)
    hold on
    for mm=1:size(S.LDBplot,1)
        LDBplotval = landboundary('read',S.LDBplot{mm,1});
        hps6(mm)=plot(LDBplotval(:,1)-S.XYoffset(1),LDBplotval(:,2)-S.XYoffset(2),S.LDBplot{mm,3},'linewidth',1.5,'MarkerSize',10,'MarkerIndices',1:10:length(x_mc0));
    end
    for sl=1:size(xp_mc,1)
        hps6(sl+mm)=plot(xp_mc{sl,:},yp_mc{sl,:},S.SLplot{sl,3},'linewidth',1.5);
        L2Name(end+1)=S.SLplot(sl,2);
        
        if S.extract_x_y
            xp{:}=xp_mc{sl,:};
            yp{:}=yp_mc{sl,:};
            out=[xp{:}',yp{:}'];
            %         namefl=cell2mat(S.SLplot(sl,2));
            save((strcat('xhrd_',num2str(S.x_hard),'_yhrd_',num2str(S.y_hard),cell2mat(S.SLplot(sl,2)),'_phi0_2',num2str(S.phiw0))),'out','-ascii');
%                     save((strcat('xhrd_','_yhrd_',cell2mat(S.SLplot(sl,2)),'_phi0_2',num2str(S.phiw0))),'out','-ascii');

        else
            out=[];
        end
    end
    hlegs = legend(hps6,L2Name','Location',S.legendlocation);
    set(hlegs,'Box','off','Color','None');
    plot(x_hard,y_hard,'k','linewidth',2);
    plot(x_mc0,y_mc0,'b--','linewidth',2);
    hold off;
    xlim(S.xlimits);
    ylim(S.ylimits);
    set(gca,'XtickLabel',num2str(get(gca,'Xtick')'/1000,'%2.1f'));
    set(gca,'YtickLabel',num2str(get(gca,'Ytick')'/1000,'%2.1f'));
    xlabel('Easting [km]');
    ylabel('Northing [km]');
    %     print((strcat('xhrd_',num2str(S.x_hard),'_yhrd_',num2str(S.y_hard),cell2mat(S.SLplot(sl,2)),'_phi0_',num2str(S.phiw0))),'-dpng', '-r300')
%     print('test1','-dpng', '-r300')
%     close(figure(12))
end
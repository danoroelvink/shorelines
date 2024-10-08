function extract_shoreline(S,STRUC,COAST,FORMAT)
% function extract_shoreline(S,STRUC,COAST,FORMAT)
% 
% For plotting the shorelines at different time steps
% export the shorelines coordinates
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

if FORMAT.printfig
    if  S.diffraction==1
%         print((strcat(S.trform,'_',num2str(2.33),'_ds0_',num2str(S.ds0),'_phi_',num2str(S.phiw0))),'-dpng', '-r300')
        print((strcat('xhrd_',num2str(STRUC.xhard),'_yhrd_',num2str(STRUC.yhard),'_phi0_2',num2str(S.phiw0))),'-dpng', '-r300')
    else  
        %print((strcat(S.trform,'_',num2str(S.b/1e4),'_ds0_',num2str(S.ds0),'_phi_',num2str(S.phiw0))),'-dpng', '-r300')
        print('-djpeg','-r300',[FORMAT.outputdir,filesep,'Final']);
    end
end

if ~isempty(FORMAT.slplot)
    if ~isempty(FORMAT.ldbplot)
        L2Name(:)=FORMAT.ldbplot(:,2);
    else
        L2Name(:)={};
        mm=0;
    end
    figure(12)
    hold on
    hps6=[];
    for mm=1:size(FORMAT.ldbplot,1)
        ldbplotval = get_landboundary(FORMAT.ldbplot{mm,1});
        hps6(mm)=plot(ldbplotval(:,1)-FORMAT.xyoffset(1),ldbplotval(:,2)-FORMAT.xyoffset(2),FORMAT.ldbplot{mm,3},'linewidth',1.5,'MarkerSize',10,'MarkerIndices',1:10:length(COAST.x_mc0));
    end
    for sl=1:size(FORMAT.xp_mc,1)
        hps6(sl+mm)=plot(FORMAT.xp_mc{sl,:},FORMAT.yp_mc{sl,:},FORMAT.slplot{sl,3},'linewidth',1.5);
        L2Name(end+1)=FORMAT.slplot(sl,2);
        
        if S.extractxy
            xp{:}=FORMAT.xp_mc{sl,:};
            yp{:}=FORMAT.yp_mc{sl,:};
            out=[xp{:}',yp{:}'];
            if S.da==1 || S.bs==1
                namefl=[cell2mat(FORMAT.slplot(sl,4)) cell2mat(FORMAT.slplot(sl,2))];
            else
                namefl=cell2mat(FORMAT.slplot(sl,2));
            end
            save(namefl,'out','-ascii');
            % save((strcat('xhrd_','_yhrd_',cell2mat(FORMAT.slplot(sl,2)),'_phi0_2',num2str(S.phiw0))),'out','-ascii');
            % save((strcat('xhrd_',num2str(STRUC.xhard),'_yhrd_',num2str(STRUC.yhard),cell2mat(FORMAT.slplot(sl,2)),'_phi0_2',num2str(S.phiw0))),'out','-ascii');
        else
            out=[];
        end
    end
    hlegs = legend(hps6,L2Name','Location',FORMAT.llocation);
    set(hlegs,'Box','off','Color','None');
    plot(STRUC.xhard,STRUC.yhard,'k','linewidth',2,'DisplayName','Structures');
    plot(COAST.x_mc0,COAST.y_mc0,'b--','linewidth',2,'DisplayName','Intial shoreline');
    hold off;
    xlim(FORMAT.xlimits);
    ylim(FORMAT.ylimits);
    set(gca,'XtickLabel',num2str(get(gca,'Xtick')'/1000,'%2.1f'));
    set(gca,'YtickLabel',num2str(get(gca,'Ytick')'/1000,'%2.1f'));
    xlabel('Easting [km]');
    ylabel('Northing [km]');
    if FORMAT.printfig
        % print((strcat('xhrd_',num2str(STRUC.xhard),'_yhrd_',num2str(STRUC.yhard),cell2mat(FORMAT.slplot(sl,2)),'_phi0_',num2str(S.phiw0))),'-dpng', '-r300')
        print('Output_without_fill','-dpng', '-r300')        
        % print('test1','-dpng', '-r300')
        % close(figure(12))
    end
end

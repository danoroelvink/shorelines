function [] = extract_berm(S,step,qwind,qwave,time,bermW,WBplot2,ds_cl,CSplot2,n_trans)
% function [] = extract_berm(S,step,qwind,qwave,time,bermW,WBplot2,ds_cl,CSplot2,n_trans)
%
% Summary of this function goes here
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

    if S.extract_berm_Plot
        save('step.txt','step','-ascii');
        save('qwind.txt','qwind','-ascii');
        save('qwave.txt','qwave','-ascii');
        save('time.txt','time','-ascii');
        saveas(figure(11),'plan.fig');
        saveas(figure(13),'trans.fig');
        if ~isempty(WBplot2)
            %[~]=plot_wb(n_trans,bermW,WBplot2);
            for i =1:n_trans
               figure(14)
               hold on
               subplot(n_trans,1,i)
               plot(WBplot2,bermW(i,:),'-b')
               datetick('x','mmm-yyyy');
               xlim([WBplot2(1) WBplot2(end)] )
               title(['Berm width variation at Transect=',sprintf('%d',i)])
               xlabel('Time');
               ylabel(' Berm Width [m]');
               print(gcf,'berm_width.jpg','-dpng','-r300')
               grid on
            end
            
            save('bermW.txt','bermW','-ascii');
            save('WBplot2.txt','WBplot2','-ascii');
        end
        if ~isempty(CSplot2)
            %[~]=plot_cl(n_trans,ds_cl,CSplot2);
            figure(14)
            hold on
            for i =1:n_trans
                subplot(n_trans,1,i)
                plot(CSplot2,ds_cl(i,:),'-b')
                datetick('x','mmm-yyyy');
                xlim([CSplot2(1) CSplot2(end)] )
                title(['Berm width variation at Transect=',sprintf('%d',i)])
                xlabel('Time');
                ylabel(' Berm Width [m]');
                grid on
            end
            print(gcf,'berm_width.jpg','-dpng','-r300')
            save('ds_cl.txt','ds_cl','-ascii');
            save('CSplot2.txt','CSplot2','-ascii');
        end
    end
end

function [FORMAT]=plot_profiles(FORMAT,TIME,DUNE,O)
% function [V,FORMAT,TIME]=plot_profiles(CHANNEL,STRUC,COAST,DUNE,WAVE,TIME,TRANSP,FORMAT,V,FNOUR)
%
% INPUT:
%   FORMAT
%      .tplot
%   O
%      .xyprofiles
%      .
%
%
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

    FORMAT.profilefighandle=22;
    if mod(TIME.it,FORMAT.plotinterval)==0 && ~isempty(O.xyprofiles)
        np=size(O.c_profile,1);
        nm=round(sqrt(np));
        nn=ceil(np./nm);
%         if isempty(FORMAT.plotprofiles)
            hf=gcf;
            figure(FORMAT.profilefighandle);
            plot_figureproperties(FORMAT.profilefighandle,900,900,32,800,50);
            if O.it_profile==0
                clf;
            end
            clear hs hp           
            for pp=1:np
                hs(pp)=subplot(nm,nn,pp);
                hp(pp,1)=plot(O.timenum_profile,O.c_profile(pp,:),'b.');
                hold on;
                if DUNE.used
                hp(pp,2)=plot(O.timenum_profile,O.d_profile(pp,:),'r.');
                end
                datetick('x','dd-mm-yyyy');
                ht=title(['profile ',num2str(pp)]);
            end
            pos0=get(gca,'Position');
            hleg=legend([hp(pp,1),hp(pp,2)],{'Coast','Dune'},'Location','SouthOutside','Orientation','Horizontal');
            set(gca,'Position',pos0);
            FORMAT.plotprofiles={hs,hp};
                
            figure(FORMAT.mainfighandle);
            %figure(hf);
%         else
%             hf=gcf;
%             figure(FORMAT.profilefighandle);
%             hs=FORMAT.plotprofiles{1};
%             hp=FORMAT.plotprofiles{2};
%             for pp=1:np
%                 set(hp(pp,1),'XData',O.timenum_profile);
%                 set(hp(pp,1),'YData',O.c_profile(pp,:));
%                 if DUNE.used
%                 set(hp(pp,2),'XData',O.timenum_profile);
%                 set(hp(pp,2),'YData',O.d_profile(pp,:));
%                 end
%             end
%             refresh(gcf);
%             refreshdata(hs(1))
%             figure(hf);
%         end
    end
end

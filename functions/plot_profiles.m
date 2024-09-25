function [FORMAT]=plot_profiles(FORMAT,TIME,DUNE,O)
% function [FORMAT]=plot_profiles(FORMAT,TIME,DUNE,O)
%
% This function plots cross-shore changes of the coastline and 
% dunefoot position over time for predefined transects along the coast. 
%
% INPUT: 
%   FORMAT
%      .plotinterval     : interval for plotting
%      .mainfighandle    : handle of the main plot figure
%      .profilefighandle : handle of the profile plot figure
%   DUNE
%      .used             : switch for using dunes (0/1)
%   TIME
%      .it               : number of timesteps since t0
%   O
%      .xyprofiles       : x and y coordinates of cross-shore transects for plotting dunes [Mx2] [m] 
%      .c_profile        : coastline position [MxNt] (m)
%      .d_profile        : dune position [MxNt] (m)
%      .it_profile       : time steps for plotting transects [1xNt] (it_profile=0 at t0)
%      .timenum_profile  : time array [days in datenum format]
% 
% OUTPUT:
%   FORMAT
%      .plotprofiles     : cell with handles of the axis (hs) and of the plotted lines (hp) as {hs,hp}
%      .profilefighandle : handle of the profile plot figure
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

    FORMAT.profilefighandle=22;
    if mod(TIME.it,FORMAT.plotinterval)==0 && ~isempty(O.xyprofiles)
        np=size(O.c_profile,1);
        nm=round(sqrt(np));
        nn=ceil(np./nm);
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

    end
end

function [S,BATHY]=update_bathy(S,BATHY,TIME,COAST)
% function [S,BATHY]=update_bathy(S,BATHY,TIME.it,COAST.x_mc,COAST.y_mc)
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

    if ~isempty(S.bathy_update)&& S.dt==(BATHY.tupdate-tnow)/365
        %% UPDATE BATHYMETRY
        S.bathyname=strcat('modified_it',num2str(TIME.it),'.dep');
        grid=wlgrid('read',S.gridfile);
        x=grid.X;
        y=grid.Y;
        S.xg=x;
        S.yg=y;
        nx=size(x,1)-1;
        ny=size(x,2)-1;
        zb=-wldep('read',S.depfile,grid);
        zb=zb(1:end-1,1:end-1);
        S.zg=zb;
        A=0.1;  %Dean profile factor 
        Dc=S.d; %closure depth

        %     figure;
        %     pcolor(S.xg,S.yg,S.zg);
        %     shading flat;
        %     colormap jet;
        %     axis equal;
        %     colorbar

        %% old function
        xg=S.xg;
        yg=S.yg;
        zg0=zb;
        %for i=1:size(xg,2);
        %    for j=1:size(xg,1);
        %        [dist(j,i),ip]=get_disttopolyline(COAST.x_mc,COAST.y_mc,xg(j,i),yg(j,i),100000); % <- 10000 in code of DR
        %    end
        %end
        zgu=zeros(size(xg));
        for i=1:size(xg,2);
            for j=1:size(xg,1);
                [dist(j,i),ip]=get_disttopolyline(COAST.x_mc,COAST.y_mc,xg(j,i),yg(j,i),100000);
               zgu(j,i)=-A*(dist(j,i)).^(2/3);
                if dist(j,i) <= S.surf_width
                   zg(j,i)=zgu(j,i);
                else
                    zg(j,i)=min(Dc,zg0(j,i));
                end
            end
        end
        IN=inpolygon(xg,yg,COAST.x_mc,COAST.y_mc);
        %zg=-dist*S.seaslope;
        %zg(zg<S.seamin)=zg0(zg<S.seamin);
        %zg=max(zg,zg0);
        %zg(IN)=min(dist(IN)*S.landslope,S.landmax);
        zg(IN)=min(A*(dist(IN)).^(2/3),S.landmax);
        zgu(IN)=min(A*(dist(IN)).^(2/3),S.landmax);
        zg(isnan(zg))=999;
        dps=-zg;
        % save(S.bathyname,'dps','-ascii')
        dps=[dps dps(:,end);dps(end,:) dps(end,end)];
        wldep('write',S.bathyname,'',dps)
        S.depfile=S.bathyname;
        % S.depsize=size(dps);
        % S.updatebathy=1;
        %save('x.dep','xg','-ascii')
        S.zg=zg;


        %% MAKE PLOT OF WAVES
        S.wavefile=strcat('Wave_table',num2str(TIME.it));
        [S]=waverefrac_implicit(S);
        close(figure(1));
        close(figure(2));

        %% UPDATE BATHY TIME AND INDEX
        BATHY.tbu=BATHY.tbu+1;
        BATHY.tupdate=BATHY.update_time(BATHY.tbu);

    end
end


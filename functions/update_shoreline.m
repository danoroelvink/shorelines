function [FORMAT,COAST]=update_shoreline(S,COAST,FORMAT)
% function [FORMAT.ifig]=update_shoreline(S)
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

    ds0=COAST.ds0;
    if ~isscalar(ds0)
        ds0=min(ds0(:,3));
    end
    
    if S.tide_interaction
        distmax=abs(S.seamin/S.seaslope);
        dntide_mc=COAST.x_mc*0;
        sedero=S.zg-S.zg0;
        yesplot=false;
        if yesplot
            figure(12);clf;
        end
        xg=S.xg;
        yg=S.yg;
        zg=S.zg;
        for i=1:size(xg,2)
            for j=1:size(xg,1)
                [dist(j,i),ip]=get_disttopolyline(COAST.x_mc,COAST.y_mc,xg(j,i),yg(j,i),10000);
                if dist(j,i)<distmax
                    dV=(sedero(j,i))*S.dsdn(j,i);
                    dntide_mc(ip)=dntide_mc(ip)+dV/ds0/COAST.h0/COAST.nt;
                    plot([COAST.x_mc(ip),xg(j,i)],[COAST.y_mc(ip),yg(j,i)],'m');
                    hold on
                end
            end
        end
        n_mc=length(find(isnan(COAST.x_mc)))+1;
        x_mc_old=COAST.x_mc;
        y_mc_old=COAST.y_mc;
        for i_mc=1:n_mc
            x=get_one_polyvar(COAST.x_mc,i_mc);
            y=get_one_polyvar(COAST.y_mc,i_mc);
            dntide=get_one_polyvar(dntide_mc,i_mc);

            n=length(x)-1;

            %% Cyclic or not ?
            cyclic = hypot(x(end)-x(1),y(end)-y(1))<ds0;
            for i=1:n
                if cyclic
                    im1=mod2(i-1,n);
                    ip1=mod2(i+1,n);
                else
                    im1=max(i-1,1);
                    ip1=min(i+1,n+1);
                end
                dn=dntide(i);
                dx(i)=-dn*(y(ip1)-y(im1))/hypot(x(ip1)-x(im1),y(ip1)-y(im1));
                dy(i)= dn*(x(ip1)-x(im1))/hypot(x(ip1)-x(im1),y(ip1)-y(im1));
            end
            for i=1:n
                x(i)=x(i)+dx(i);
                y(i)=y(i)+dy(i);
            end
            if cyclic
                x(n+1)=x(1);
                y(n+1)=y(1);
            end

            [COAST.x_mc,COAST.y_mc]=insert_section(x,y,COAST.x_mc,COAST.y_mc,i_mc);
        end

        COAST.dntide=dntide_mc;
        pcolor(xg,yg,sedero);
        cmax=max(max(abs(sedero)));
        caxis([-cmax cmax])
        colorbar;
        shading flat;
        hold on
        plot(x_mc_old,y_mc_old,COAST.x_mc,COAST.y_mc,'linewidth',2);
        axis equal
        if isfield(FORMAT,'ifig')
            FORMAT.ifig=FORMAT.ifig+1;
        else
            FORMAT.ifig=1000;
        end
        % fname=[S.outputdir,filesep,'sedero',num2str(S.ifig),'.jpg'];
        % print('-djpeg',fname)

        %%pause
    end
end

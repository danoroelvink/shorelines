function [CLplot,iint,CLplot2,ds_cl]=plot_coastline_transect(S,x_mc,y_mc,x_mc0,y_mc0,n_trans,CLplot,CLplot2,tnow,iint,x_trans,y_trans,ds_cl)
% function [CLplot,iint,CLplot2,ds_cl]=plot_coastline_transect(S,x_mc,y_mc,x_mc0,y_mc0,n_trans,CLplot,CLplot2,tnow,iint,x_trans,y_trans,ds_cl)
% 
% a function to plot transects of cross-shore distances for both coastline and dune foot to a reference line...
% at specific times and transecs
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


    if S.dt==(CLplot-tnow)/365
        for i =1:n_trans
            xs_trans=[x_trans(2*i-1) x_trans((2*i-1)+1)];
            ys_trans=[y_trans(2*i-1) y_trans((2*i-1)+1)];
            [xx1,yy1]=get_intersections(xs_trans,ys_trans,x_mc0,y_mc0);
            [xx2,yy2]=get_intersections(xs_trans,ys_trans,x_mc,y_mc);
            ds=hypot(XX1(1)-XX2(1),yy1(1)-yy2(1));
            ds_cl(i,iint)=ds;
        end
        CLplot2(iint)=CLplot;
        CLplot_t=datetime(CLplot,'ConvertFrom','datenum','Format','yyyy-MM-dd');
        CLplot=CLplot_t+calmonths(S.CLplot_int);
        CLplot=datenum(CLplot);
        iint=iint+1;
    end
end

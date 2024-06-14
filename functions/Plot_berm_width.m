function [BWplot, bermW,int,BWplot2]=plot_berm_width(S,x_mc,y_mc,x_dune,y_dune,xl,yl,n_trans,BWplot,BWplot2,bermW,tnow,int,x_trans,y_trans)
% function [BWplot, bermW,int,BWplot2]=plot_berm_width(S,x_mc,y_mc,x_dune,y_dune,xl,yl,n_trans,BWplot,BWplot2,bermW,tnow,int,x_trans,y_trans)
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


    if S.dt==(BWplot-tnow)/365
        for i =1:n_trans
            xs_trans=[x_trans(2*i-1) x_trans((2*i-1)+1)];
            ys_trans=[y_trans(2*i-1) y_trans((2*i-1)+1)];
            [xx1,yy1]=get_intersections(xs_trans,ys_trans,xl,yl);
            [xx2,yy2]=get_intersections(xs_trans,ys_trans,x_mc,y_mc);
            [xx3,yy3]=get_intersections(xs_trans,ys_trans,x_dune,y_dune);
            dis_cl=hypot(xx1(1)-xx2(1),yy1(1)-yy2(1));
            dis_dn=hypot(xx1(1)-xx3(1),yy1(1)-yy3(1));
            bermw=dis_cl-dis_dn;
            bermW(i,int)=bermw;
        end
        BWplot2(int)=BWplot;
        BWplot_t=datetime(BWplot,'ConvertFrom','datenum','Format','yyyy-MM-dd');
        BWplot=BWplot_t+calmonths(S.bermw_plot_int);
        BWplot=datenum(BWplot);
        int=int+1;
    end
end

function [xl,yl,phirf,vi,vii,xmax,ymax,xmin,ymin,nmax,iwtw,plot_time,tsl,tplot,xp_mc,yp_mc,CLplot,CLplot2,iint,BWplot,BWplot2,innt,ds_cl,qwave,qwind,time,step,int,bermW,x_trans,y_trans,n_trans] = initialize_plot_variables(S,x_mc,y_mc,x_hard,y_hard,x_dune,y_dune,timenum0)
% function [xl,yl,phirf,vi,vii,xmax,ymax,xmin,ymin,nmax,iwtw,plot_time,tsl,tplot,xp_mc,yp_mc,CLplot,CLplot2,iint,BWplot,BWplot2,innt,ds_cl,qwave,qwind,time,step,int,bermW,x_trans,y_trans,n_trans] = initialize_plot_variables(S,x_mc,y_mc,x_hard,y_hard,x_dune,y_dune,timenum0)
%
% UNTITLED5 Summary of this function goes here
% Detailed explanation goes here
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

    %  Create reference line for cross-shore plotting
    [xl,yl,phirf]= dune_refline(S,x_dune,y_dune);
    
    % for making video
    vi = struct('cdata', cell(1,1), 'colormap', cell(1,1));
    vii=0;
    
    %for plotting fill
    xmax=0;ymax=0;xmin=0;ymin=0;nmax=0;
    iwtw=0;  % for using interpolated wave table (warning)
    
    % For extracting specific shorelines
    if ~isempty(S.SLplot)
        plot_time(:)=datenum(S.SLplot(:,1),'yyyy-mm-dd'); % HH:MM:SS
        tsl=1;
        tplot=plot_time(tsl);
        plot_time(end+1)=0;
    else
        tplot=[];
        tsl=[];
        plot_time=[];
    end
    xp_mc{:,:}={};
    yp_mc{:,:}={};
    if S.dune && S.qplot                                                         % for extracting wave and wind transport rates at a certain transect(s)
        innt=1;
    else
        innt=[];
    end
    if S.dune && S.bermw_plot 
        BWplot=timenum0;                                                       % for extracting beach berm width at a certain transect(s)
        BWplot2=timenum0;
        int=1;
    else
        BWplot=[];
        BWplot2=[];
        int=[];
    end
    if S.CLplot
        CLplot=timenum0;                                                       % for extracting coastline position relative to the initial coastline at a certain transect(s)
        CLplot2=timenum0;
        iint=1;
    else
        CLplot=[];
        CLplot2=[];    
        iint=1;
    end
    ds_cl=[];
    qwave=[];
    qwind=[];
    time=[];
    step=[];
    bermW=[];
    
    %% prepare transections from cross-shore plotting
    if S.CLplot ||S.bermw_plot ||S.qplot 
        [x_trans,y_trans,n_trans]=prepare_transects(S,x_mc,y_mc,x_hard,y_hard,x_dune,y_dune,xl,yl);
    else
        x_trans=[];
        y_trans=[];
        n_trans=[];
    end
end

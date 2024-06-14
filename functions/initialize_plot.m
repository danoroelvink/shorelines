function [FORMAT] = initialize_plot(S,COAST)
% function [FORMAT.plot_time,FORMAT.tsl,FORMAT.tplot,FORMAT.xp_mc,FORMAT.yp_mc] = initialize_plot_variables(S)
% function [xl,yl,phirf,vi,vii,FORMAT.xmax,FORMAT.ymax,FORMAT.xmin,FORMAT.ymin,FORMAT.nmax,iwtw,FORMAT.plot_time,FORMAT.tsl,FORMAT.tplot,FORMAT.xp_mc,FORMAT.yp_mc,CLplot,CLplot2,iint,BWplot,BWplot2,innt,ds_cl,qwave,qwind,time,step,int,bermW,x_trans,y_trans,n_trans] = initialize_plot_variables(S,x_mc,y_mc,x_hard,y_hard,x_dune,y_dune,timenum0)
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

    fprintf('  Initialize plot \n');
       
    % Set x and y limits
    if isempty(S.xlimits)
        xlims=[min(COAST.x_mc),max(COAST.x_mc)];
        S.xlimits = [min(xlims)-0.05*diff(xlims),max(xlims)+0.05*diff(xlims)] - S.XYoffset(1);
    else
        S.xlimits = S.xlimits - S.XYoffset(1);
    end
    if isempty(S.ylimits)
        ylims=[min(COAST.y_mc),max(COAST.y_mc)];
        S.ylimits = [min(ylims)-0.05*diff(ylims),max(ylims)+0.05*diff(ylims)] - S.XYoffset(2);
    else
        S.ylimits = S.ylimits - S.XYoffset(2);
    end
    S.xlimits(2)=max(S.xlimits(1)+500,S.xlimits(2)); 
    S.ylimits(2)=max(S.ylimits(1)+500,S.ylimits(2)); 
    S.xlimits=sort(S.xlimits); 
    S.ylimits=sort(S.ylimits); 
    FORMAT.xlimits=S.xlimits; 
    FORMAT.ylimits=S.ylimits; 
    FORMAT.XYoffset=S.XYoffset; 

    % Format figures
    if S.plotvisible==0
        FORMAT.plotvisible='off';
    else
        FORMAT.plotvisible='on';
    end
    
    FORMAT.mainfighandle = figure(11);clf;
    set(FORMAT.mainfighandle,'Visible',FORMAT.plotvisible);
    FORMAT.xmax=0;FORMAT.ymax=0;FORMAT.xmin=0;FORMAT.ymin=0;FORMAT.nmax=0;
    % set properties of figure
    if diff(S.xlimits)/diff(S.ylimits) > 1850/850
        sx=1850;
        sy=1850*diff(S.ylimits)/diff(S.xlimits);
    else
        sx=850*diff(S.xlimits)/diff(S.ylimits);
        sy=850;
    end
    plot_figureproperties(FORMAT.mainfighandle,sx,sy,32);
    xlim(S.xlimits);ylim(S.ylimits); 
    
    try
        FORMAT.IMGfilename=S.IMGfilename;
        FORMAT.WORLDfilename=S.WORLDfilename;
    end
    
    if S.yesplot
        figure(11);clf;
        plot([S.xlimits(1) S.xlimits(2) S.xlimits(2) S.xlimits(1) S.xlimits(1)], ...
            [S.ylimits(1) S.ylimits(1) S.ylimits(2) S.ylimits(2) S.ylimits(1)],'k:');
        axis equal;
        hold on;
        xlabel('Easting [m]');
        ylabel('Northing [m]');
        plot(COAST.x_mc0,COAST.y_mc0,'k','linewidth',2);
        axis equal
        if ~isempty(S.xlimits)
            xlim(S.xlimits);
        end
        if ~isempty(S.ylimits)
            ylim(S.ylimits);
        end
    end

    FORMAT.xp_mc={};
    FORMAT.yp_mc={};
    
    % empty land polygon 
    FORMAT.xb=[];
    FORMAT.yb=[];
    
    % Other plot formatting
    if ~isempty(S.SLplot)  % For extracting specific shorelines
        FORMAT.plot_time(:)=datenum(S.SLplot(:,1),'yyyy-mm-dd');
        FORMAT.tsl=1;
        FORMAT.tplot=FORMAT.plot_time(FORMAT.tsl);
        FORMAT.plot_time(end+1)=0;
    else
        FORMAT.tplot=[];
        FORMAT.tsl=[];
        FORMAT.plot_time=[];
    end
    
    % create output directory
    if ~exist(fullfile(pwd,S.outputdir),'dir')
        mkdir(fullfile(pwd,S.outputdir));
    end
    
    FORMAT.plotinterval=S.plotinterval;
    FORMAT.usefill=S.usefill;
    FORMAT.usefillpoints=S.usefillpoints;
    FORMAT.LDBplot=S.LDBplot;
    FORMAT.llocation=S.llocation;
    FORMAT.SLplot=S.SLplot;
    FORMAT.video=S.video;
    FORMAT.outputdir=S.outputdir;
    FORMAT.fignryear=S.fignryear;
    FORMAT.ld=S.ld;
    FORMAT.XYwave=S.XYwave;
    FORMAT.plotQS=S.plotQS;
    FORMAT.plotHS=S.plotHS;
    FORMAT.plotDIR=S.plotDIR;
    FORMAT.plotUPW=S.plotUPW;
    FORMAT.print_fig=S.print_fig;
    FORMAT.fastplot=S.fastplot;
    FORMAT.plotprofiles=[];
    FORMAT.xyprofiles=S.xyprofiles;
    FORMAT.plotcoast=0;
    %% write logfile
    % struct2log(FORMAT,'FORMAT','a');

end
 
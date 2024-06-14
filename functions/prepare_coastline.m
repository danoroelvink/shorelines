function [COAST]=prepare_coastline(S)
% function [x_mc,y_mc,x_mc0,y_mc0,S]=prepare_coastline(S)
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

    fprintf('  Prepare coastline \n');
    COAST=struct;
    COAST.LDBcoastline=S.LDBcoastline;
    COAST.h0=[];
    COAST.h0input=S.d;
    COAST.ds0=S.ds0;
    COAST.xlimits=S.xlimits;
    COAST.ylimits=S.ylimits;
    COAST.XYoffset=S.XYoffset;
    COAST.twopoints=S.twopoints;
    COAST.smoothfac=S.smoothfac;
    COAST.smoothrefrac=S.smoothrefrac;
    COAST.PHIf0=S.phif;
    COAST.tanbeta=S.tanbeta;
    COAST.griddingmethod=S.griddingmethod;
    COAST.maxangle=S.maxangle;
    
    % boundary conditions
    % {'Closed',e.g. 0 or 9000 m3/yr}, {'Neumann',dummy},{'Fixed',dummy},{'Angleconstant',empty to use at t0 or specified value e.g. 321°N};
    COAST.PHIc0bnd=[nan,nan];
    if ischar(S.boundary_condition_start)
        COAST.boundary_condition_start={S.boundary_condition_start,nan};                     % boundary condition 'Closed', 'Neumann','Fixed','Angleconstant'
    else
        COAST.boundary_condition_start=S.boundary_condition_start;
    end
    if ischar(S.boundary_condition_end) 
        COAST.boundary_condition_end={S.boundary_condition_end,nan};                       % boundary condition 'Closed', 'Neumann','Fixed','Angleconstant'
    else
        COAST.boundary_condition_end=S.boundary_condition_end;
    end    
    
    if ~isempty(S.LDBcoastline) && ~isfield(S,'x_mc')
        if ischar(S.LDBcoastline)
            xy_mc=load(S.LDBcoastline);
        else
            xy_mc=S.LDBcoastline;
        end
        if isempty(S.XYoffset)
            COAST.XYoffset = [floor(min(xy_mc(:,1))/1000)*1000 , floor(min(xy_mc(:,2))/1000)*1000];
        end
        COAST.x_mc=xy_mc(:,1)' - S.XYoffset(1);   %COAST.x_mc=xy_mc(end:-1:1,1)';   % SHIFT COASLTINE
        COAST.y_mc=xy_mc(:,2)' - S.XYoffset(2);   %COAST.y_mc=xy_mc(end:-1:1,2)';   % SHIFT COASLTINE
        COAST.x_mc0=COAST.x_mc;
        COAST.y_mc0=COAST.y_mc;

    elseif ~isfield(S,'x_mc')
        figure(11);clf;
        %plot_figureproperties(gcf,800,950,32);
        %plot_figureproperties(gcf,920,1120,32,80,0);
        %plot_figureproperties(gcf,1855,1120,32,65,0);
        xl=S.xlimits;yl=S.ylimits;
        plot([S.xlimits(1) S.xlimits(2) S.xlimits(2) S.xlimits(1) S.xlimits(1)], ...
            [S.ylimits(1) S.ylimits(1) S.ylimits(2) S.ylimits(2) S.ylimits(1)],'k:');
        axis equal;
        %xl=xlim;yl=ylim;
        xlabel('Easting [m]');
        ylabel('Northing [m]');
        htxt=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add coastline (LMB); Next segment (RMB); Exit (q)');set(htxt,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.1 0.6]);
        [COAST.x_mc,COAST.y_mc]=select_multi_polygon('r');
        set(htxt,'Visible','off');
        COAST.x_mc0=COAST.x_mc;
        COAST.y_mc0=COAST.y_mc;

    else
        COAST.x_mc=S.x_mc;
        COAST.y_mc=S.y_mc;
        try
            COAST.x_mc0=S.x_mc0;
            COAST.y_mc0=S.y_mc0;
        catch 
            COAST.x_mc0=S.x_mc;
            COAST.y_mc0=S.y_mc;
        end
        figure(11);clf;
    end
    
    % tidy up the coastline without nans at the start and end, and remove double nan's
    idnotnan=find(~isnan(COAST.x_mc));
    COAST.x_mc=COAST.x_mc(idnotnan(1):idnotnan(end));
    COAST.y_mc=COAST.y_mc(idnotnan(1):idnotnan(end));
    idnan=find(isnan(COAST.x_mc));
    iduse=setdiff([1:length(COAST.x_mc)],idnan(diff(idnan)==1));
    COAST.x_mc=COAST.x_mc(iduse);
    COAST.y_mc=COAST.y_mc(iduse);

    %% Invert order of coastline points if defined in GIS conventions (e.g. ccw)
    if S.gis_convention
        COAST.x_mc = COAST.x_mc(end:-1:1);
        COAST.y_mc = COAST.y_mc(end:-1:1);
    end

    % determine number of coastline segments
    nans=find(isnan(COAST.x_mc));
    COAST.n_mc=length(nans)+1;
    
    % coast angles (PHIcxy_mc0) at coastline points only for plotting DUNE
    xq=(COAST.x_mc0(1:end-1)+COAST.x_mc0(2:end))/2;
    yq=(COAST.y_mc0(1:end-1)+COAST.y_mc0(2:end))/2;
    xq=[1.5*COAST.x_mc0(1)-0.5*COAST.x_mc0(2),xq];
    yq=[1.5*COAST.y_mc0(1)-0.5*COAST.y_mc0(2),yq];
    xq=[xq,COAST.x_mc0(end)];
    yq=[yq,COAST.y_mc0(end)];
    COAST.PHIcxy_mc0 = mod(360.0-atan2d(diff(yq),diff(xq)),360);
    
    %% write logfile
    % struct2log(COAST,'COAST','a');

end 

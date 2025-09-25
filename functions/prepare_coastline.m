function [COAST]=prepare_coastline(S)
% function [COAST]=prepare_coastline(S)
%
% This function reads all the input variables and initializes the internal data-structure COAST.
% 
% INPUT:
%        .ldbcoastline (option 1) : filename with coastline data [Nx2]
%        .x_mc         (option 2) : x-coordinates of initial coastline [m]
%        .y_mc         (option 2) : y-coordinates of initial coastline [m]
%        <other, see output>
% 
% OUTPUT:
%    COAST
%        .x_mc                    : x-coordinates of coastal points for all coastal elements (m) [1xN]
%        .y_mc                    : y-coordinates of coastal points for all coastal elements (m) [1xN]
%        .n_mc                    : number of output points
%        .x_mc0                   : x-coordinates of initial coastline at t0 [m]
%        .y_mc0                   : y-coordinates of initial coastline at t0 [m]
%        .h0                      : active height along the coast (as used) [m]
%        .h0input                 : active height input [m] (e.g. value or [x1,y1,h1; etc] )
%        .ds0                     : grid cell size [m]
%        .xlimits                 : xlimits used for plotting [m]
%        .ylimits                 : ylimits used for plotting [m]
%        .xyoffset                : xy-offset used for plotting (1x2) [m]
%        .twopoints               : switch for 'S.twopoints approach' which determines the type of response during high-angle wave events (default=1)
%        .smoothfac               : smoothing factor used to re-arrange grid every timestep (only for griddingmethod==1), with reasonable values between 0 and 0.1.
%        .smoothrefrac            : smoothing fraction of the coastline orientation (PHIcs) which is used only for the refraction of waves from nearshore location (TDP) to point of breaking (BR), and not for the transport computation. Values between 0 and 1.
%        .PHIf0                   : Orientation of the foreshore [°N]
%        .tanbeta                 : mean bed slope [ratio 1/slope] 
%        .griddingmethod          : method for regenerating the grid (1: only splitting and merging cells if they are too small, 2: uniform grid regeneration if criteria for gridsize or exceeded)
%        .maxangle                : maximum coastline re-orientation between individual grid cells (affecting spit width and stabilizing small scale features in case of dense grids)
%        .PHIc0bnd                : coastline angle at boundary points at t0
%        .boundaryconditionstart  : boundary condition at the start of the coastal element (e.g. 'Closed', 'Neumann','Fixed' or {'Angleconstant',145} )
%        .boundaryconditionend    : boundary condition at the end of the coastal element (e.g. 'Closed', 'Neumann','Fixed' or {'Angleconstant',145} )
%        .PHIcxy_mc0              : orientation of the coastline at t0 [°N]
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

    fprintf('  Prepare coastline \n');
    COAST=struct;
    COAST.ldbcoastline=S.ldbcoastline;
    COAST.h0=[];
    COAST.h0input=S.d;
    COAST.ds0=S.ds0;
    COAST.xlimits=S.xlimits;
    COAST.ylimits=S.ylimits;
    COAST.xyoffset=S.xyoffset;
    COAST.twopoints=S.twopoints;
    COAST.smoothfac=S.smoothfac;
    COAST.smoothrefrac=S.smoothrefrac;
    COAST.PHIf0=S.phif;
    COAST.tanbeta=S.tanbeta;
    COAST.griddingmethod=S.griddingmethod;
    COAST.maxangle=S.maxangle;
    COAST.preserveorientation=S.preserveorientation;
    COAST.mergegrid=0;
    
    % boundary conditions
    % {'Closed',e.g. 0 or 9000 m3/yr}, {'Neumann',dummy},{'Fixed',dummy},{'Angleconstant',empty to use at t0 or specified value e.g. 321°N};
    COAST.PHIc0bnd=[nan,nan];
    if ischar(S.boundaryconditionstart)
        COAST.boundaryconditionstart={S.boundaryconditionstart,nan};                     % boundary condition 'Closed', 'Neumann','Fixed','Angleconstant'
    else
        COAST.boundaryconditionstart=S.boundaryconditionstart;
    end
    if ischar(S.boundaryconditionend) 
        COAST.boundaryconditionend={S.boundaryconditionend,nan};                       % boundary condition 'Closed', 'Neumann','Fixed','Angleconstant'
    else
        COAST.boundaryconditionend=S.boundaryconditionend;
    end    
    
    if ~isempty(S.ldbcoastline)
        % Read coastline from S.ldbcoastline
        if ischar(S.ldbcoastline)
            xy_mc=load(S.ldbcoastline);
        else
            xy_mc=S.ldbcoastline;
        end
        if isempty(S.xyoffset)
            COAST.xyoffset = [floor(min(xy_mc(:,1))/1000)*1000 , floor(min(xy_mc(:,2))/1000)*1000];
        end
        COAST.x_mc=xy_mc(:,1)' - S.xyoffset(1);   %COAST.x_mc=xy_mc(end:-1:1,1)';   % SHIFT COASLTINE
        COAST.y_mc=xy_mc(:,2)' - S.xyoffset(2);   %COAST.y_mc=xy_mc(end:-1:1,2)';   % SHIFT COASLTINE
        COAST.x_mc0=COAST.x_mc;
        COAST.y_mc0=COAST.y_mc;

    elseif ~isempty(S.xmc)
        % Read coastline from S.xmc and S.ymc
        COAST.x_mc=S.xmc;
        COAST.y_mc=S.ymc;
        try
            COAST.x_mc0=S.xmc0;
            COAST.y_mc0=S.ymc0;
        catch 
            COAST.x_mc0=S.xmc;
            COAST.y_mc0=S.ymc;
        end
        figure(11);clf;

    else
        % Interactive mode
        figure(11);clf;
        xl=S.xlimits;yl=S.ylimits;
        plot([S.xlimits(1) S.xlimits(2) S.xlimits(2) S.xlimits(1) S.xlimits(1)], ...
            [S.ylimits(1) S.ylimits(1) S.ylimits(2) S.ylimits(2) S.ylimits(1)],'k:');
        axis equal;
        xlabel('Easting [m]');
        ylabel('Northing [m]');
        htxt=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add coastline (LMB); Next segment (RMB); Exit (q)');set(htxt,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.1 0.6]);
        [COAST.x_mc,COAST.y_mc]=select_multi_polygon('r');
        set(htxt,'Visible','off');
        COAST.x_mc0=COAST.x_mc;
        COAST.y_mc0=COAST.y_mc;
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
    if S.gisconvention
        COAST.x_mc = COAST.x_mc(end:-1:1);
        COAST.y_mc = COAST.y_mc(end:-1:1);
    end

    % determine number of coastline segments
    nans=find(isnan(COAST.x_mc));
    COAST.n_mc=length(nans)+1;
    COAST.n_mc0=COAST.n_mc;
    
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

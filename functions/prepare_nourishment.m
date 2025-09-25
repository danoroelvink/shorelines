function [NOUR] = prepare_nourishment(S,COAST,STRUC)
% function [NOUR] = prepare_nourishment(S,COAST,STRUC)
%
% The nourishment data is initialized, and the data-structure NOUR filled with variables.
%
% INPUT: 
%    S
%       .nourish         : switch for using nourishments (0/1)
%       .growth          : factor for scaling the nourishment efficiency (from 0 to 1, default is 1)
%       <option 1>
%       .norfile         : file with nourishment information with .NOR extension. Containing a row for each nourishment with [xstart, ystart,  xend,  yend,  tstart(yyyymmdd), tend(yyyymmdd), totalvolume]
%       <option 2>
%       .ldbnourish      : LDB with nourishment locations ([Nx2] ASCII FILE WITHOUT HEADER) (i.e. polygon around relevant grid cells) <- leave empty to use interactive mode!
%       .nourratefile    : nourishment rates placed in order for each nourishment polygons (not used if .NOR file is available
%       .nourstartfile   : nourishment start dates placed in order for each  nourishment polygons (not used if .NOR file is available
%       .nourendfile     : nourishment end dates placed in order for each  nourishment polygons (not used if .NOR file is available
%       .nourrate        : rate of nourishing (not used if .NOR file is available)
%       .nourmethod      : a method can be chosen which uses only the begin and end point of the nourishment to identify where the nourishment needs to take place ('default') or a method that can distribute the sediment over more than 2 elements ('complex').
%       .reftime         : reference time (i.e. 'yyyy-mm-dd')
%       .endofsimulation : end time (i.e. 'yyyy-mm-dd')
%    COAST
%       .x_mc            : x-coordinates of coastline points (only needed for plotting / interactive mode)
%       .y_mc            : y-coordinates of coastline points (only needed for plotting / interactive mode)
%    STRUC
%       .xhard           : x-coordinates of hard structures (only needed for plotting / interactive mode)
%       .yhard           : y-coordinates of hard structures (only needed for plotting / interactive mode)
%
%
% OUTPUT: 
%    NOUR
%       .nourish         : switch for using nourishments (0, 1 or 2) -> 2 is common method used
%       .growth	         : factor for scaling the nourishment efficiency (from 0 to 1, default is 1)
%       .xnour           : x-coordinates of nourishments [Nx1]
%       .ynour           : y-coordinates of nourishments [Nx1]
%       .nnour           : number of nourishments
%       .tstart          : start time of nourishments [days in datenum format]
%       .tend            : end time of nourishments [days in datenum format]
%       .rate            : nourishment rate [m3/yr]
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

    fprintf('  Prepare nourishments \n');
    NOUR=struct;
    NOUR.xnour=[];
    NOUR.ynour=[];
    NOUR.nnour=[];
    NOUR.tstart=[];
    NOUR.tend=[];
    NOUR.rate=[];
    NOUR.nourish=S.nourish;
    NOUR.growth=S.growth;
    NOUR.method=S.nourmethod;
    
    % in case keyword 'S.norfile' is used instead of 'ldbnourish'.
    if ~isempty(findstr(lower(S.norfile),'.nor'))
        S.ldbnourish=S.norfile;
        NOUR.nourish=1;
    elseif isfield(S,'ldbnour') & isempty(S.ldbnourish)
        S.ldbnourish=S.ldbnour;
        NOUR.nourish=1;
    elseif ~isempty(findstr(lower(S.ldbnourish),'.nor'))
        NOUR.nourish=1;
    end
    
    if ~isempty(NOUR.nourish) && NOUR.nourish ~= 0
        
        if ~isempty(findstr(lower(S.ldbnourish),'.nor'))
            NOUR.nourish=2;
            readmatfile=1;
            % LBDnourish contains [xstart, ystart,  xend,  yend,  tstart(yyyymmdd), tend(yyyymmdd), totalvolume]
            nor=load(S.ldbnourish);
            NOUR.nnour=size(nor,1);
            NOUR.xnour=[nor(:,1),nor(:,3)];
            NOUR.ynour=[nor(:,2),nor(:,4)];
            NOUR.tstart=datenum(num2str(nor(:,5)),'yyyymmdd');
            NOUR.tend=datenum(num2str(nor(:,6)),'yyyymmdd');
            NOUR.rate=(nor(:,7)./(NOUR.tend-NOUR.tstart))*365; % m3/year

        elseif ~isempty(S.ldbnourish) && NOUR.nourish==3
            NOUR.nourish=2;
            readmatfile=1;
            % LBDnourish contains [xstart, ystart,  xend,  yend,  tstart(yyyymmdd), tend(yyyymmdd), totalvolume]
            nor=load(S.ldbnourish);
            nor = nor.nor;
            NOUR.nnour=size(nor,1);
            NOUR.xnour=[nor(:,1),nor(:,3)];
            NOUR.ynour=[nor(:,2),nor(:,4)];
            NOUR.tstart=datenum(num2str(nor(:,5)),'yyyymmdd');
            NOUR.tend=datenum(num2str(nor(:,6)),'yyyymmdd');
            NOUR.rate=(nor(:,7)./(NOUR.tend-NOUR.tstart))*365; % m3/year

        elseif ~isempty(S.ldbnourish)
            NOUR.nourish=1;
            readmatfile=0;
            xynour=load(S.ldbnourish);
            NOUR.xnour=xynour(:,1)'-S.xyoffset(1);
            NOUR.ynour=xynour(:,2)'-S.xyoffset(2);
            matname=S.ldbnourish;
            matname(end-2:end)='mat';
            if exist(matname)==2
                load(matname);
                NOUR.nnour=length(nourishments);
                for inour=1:NOUR.nnour
                    NOUR.tstart(inour) = nourishments(inour).tstart;
                    NOUR.tend(inour)   = nourishments(inour).tend;
                    NOUR.rate(inour) = nourishments(inour).rate;
                end
                readmatfile=1;
            end

        else
            NOUR.nourish=1;
            readmatfile=0;

            if length(STRUC.xhard)>1
                figure(11);
                plot(COAST.x_mc,COAST.y_mc,'-*k',STRUC.xhard,STRUC.yhard,'-r');
                axis equal
            else
                figure(11);
                plot(COAST.x_mc,COAST.y_mc,'-*k');
                axis equal
            end
            xl=xlim;yl=ylim;
            htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add nourishment area(LMB); Next nourishment area (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
            [NOUR.xnour,NOUR.ynour]=select_multi_polygon('k');
            set(htxt2,'Visible','off');
        end
        if readmatfile==0
            [ ~,~,NOUR.nnour,~,~ ] = get_one_polygon( NOUR.xnour,NOUR.ynour,1 );

            if ~isempty(S.nourratefile)             % when using arrays of dates and nourishment rates
                NOUR.tstart=datenum(importdata(S.nourstartfile));
                NOUR.tend=datenum(importdata(S.nourendfile));
                NOUR.rate=load(S.nourratefile);
            else
                NOUR.tstart(1:n)=datenum(S.reftime);
                NOUR.tend(1:n)=datenum(S.endofsimulation);
                NOUR.rate(1:NOUR.nnour)=S.nourrate;
            end
        end
        %% write logfile
        % struct2log(NOUR,'NOUR','a');
    end

end

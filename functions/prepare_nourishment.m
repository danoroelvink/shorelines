function [NOUR] = prepare_nourishment(S,COAST,STRUC)
% function [NOUR] = prepare_nourishment(S,COAST,STRUC)
%
% UNTITLED4 Summary of this function goes here
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

    fprintf('  Prepare nourishments \n');
    NOUR=struct;
    NOUR.x_nour=[];
    NOUR.y_nour=[];
    NOUR.n_nour=[];
    NOUR.tstart=[];
    NOUR.tend=[];
    NOUR.rate_m3_per_yr=[];
    NOUR.nourish=S.nourish;
    NOUR.growth=S.growth;
    
    % in case keyword 'S.norfile' is used instead of 'LDBnourish'.
    if ~isempty(findstr(lower(S.norfile),'.nor'))
        S.LDBnourish=S.norfile;
        NOUR.nourish=1;
    end
    
    if ~isempty(NOUR.nourish) && NOUR.nourish ~= 0
        
        if ~isempty(findstr(lower(S.LDBnourish),'.nor'))
            NOUR.nourish=2;
            readmatfile=1;
            % LBDnourish contains [xstart, ystart,  xend,  yend,  tstart(yyyymmdd), tend(yyyymmdd), totalvolume]
            nor=load(S.LDBnourish);
            NOUR.n_nour=size(nor,1);
            NOUR.x_nour=[nor(:,1),nor(:,3)];
            NOUR.y_nour=[nor(:,2),nor(:,4)];
            NOUR.tstart=datenum(num2str(nor(:,5)),'yyyymmdd');
            NOUR.tend=datenum(num2str(nor(:,6)),'yyyymmdd');
            NOUR.rate_m3_per_yr=(nor(:,7)./(NOUR.tend-NOUR.tstart))*365; % m3/year

        elseif ~isempty(S.LDBnourish) && NOUR.nourish==3
            NOUR.nourish=2;
            readmatfile=1;
            % LBDnourish contains [xstart, ystart,  xend,  yend,  tstart(yyyymmdd), tend(yyyymmdd), totalvolume]
            nor=load(S.LDBnourish);
            nor = nor.nor;
            NOUR.n_nour=size(nor,1);
            NOUR.x_nour=[nor(:,1),nor(:,3)];
            NOUR.y_nour=[nor(:,2),nor(:,4)];
            NOUR.tstart=datenum(num2str(nor(:,5)),'yyyymmdd');
            NOUR.tend=datenum(num2str(nor(:,6)),'yyyymmdd');
            NOUR.rate_m3_per_yr=(nor(:,7)./(NOUR.tend-NOUR.tstart))*365; % m3/year

        elseif ~isempty(S.LDBnourish)
            NOUR.nourish=1;
            readmatfile=0;
            xy_nour=load(S.LDBnourish);
            NOUR.x_nour=xy_nour(:,1)'-S.XYoffset(1);
            NOUR.y_nour=xy_nour(:,2)'-S.XYoffset(2);
            matname=S.LDBnourish;
            matname(end-2:end)='mat';
            if exist(matname)==2
                load(matname);
                NOUR.n_nour=length(nourishments);
                for inour=1:NOUR.n_nour
                    NOUR.tstart(inour) = nourishments(inour).tstart;
                    NOUR.tend(inour)   = nourishments(inour).tend;
                    NOUR.rate_m3_per_yr(inour) = nourishments(inour).rate;
                end
                readmatfile=1;
            end   

        else
            NOUR.nourish=1;
            readmatfile=0;

            if length(STRUC.x_hard)>1
                figure(11);
                plot(COAST.x_mc,COAST.y_mc,'-*k',STRUC.x_hard,STRUC.y_hard,'-r');
                axis equal
            else
                figure(11);
                plot(COAST.x_mc,COAST.y_mc,'-*k');
                axis equal
            end
            xl=xlim;yl=ylim;
            htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add nourishment area(LMB); Next nourishment area (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
            [NOUR.x_nour,NOUR.y_nour]=select_multi_polygon('k');
            set(htxt2,'Visible','off');
        end
        if readmatfile==0
            [ ~,~,NOUR.n_nour,~,~ ] = get_one_polygon( NOUR.x_nour,NOUR.y_nour,1 );

            if ~isempty(S.nourratefile)             % when using arrays of dates and nourishment rates
                NOUR.tstart=datenum(importdata(S.nourstartfile));
                NOUR.tend=datenum(importdata(S.nourendfile));
                NOUR.rate_m3_per_yr=load(S.nourratefile);
            else
                NOUR.tstart(1:n)=datenum(S.reftime);
                NOUR.tend(1:n)=datenum(S.endofsimulation);
                NOUR.rate_m3_per_yr(1:NOUR.n_nour)=S.nourrate;
            end
        end
        %% write logfile
        % struct2log(NOUR,'NOUR','a');
    end

end

function [COAST,merged] = merge_coastlines(COAST,i_mc)
% function [COAST,merged] = merge_coastlines(COAST,i_mc)
% 
% Merges and splits individual coastal elements that cross themselves.
% For example, a spit tip reaching the beach. 
%
% INPUT:
%    COAST
%         .x_mc  : x-coordinates of the coastline elements 
%         .y_mc  : y-coordinates of the coastline elements 
%         .ds0   : grid cell size [m]
% 
% OUTPUT:
%    COAST
%         .x_mc  : x-coordinates of the coastline elements (after splitting/merging)
%         .y_mc  : y-coordinates of the coastline elements (after splitting/merging)
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

    %eps=.1;
    eps=1;
    s(1)=0;
    yesplot=0;
    [ COAST.x,COAST.y,COAST.n_mc ] = get_one_polygon( COAST.x_mc,COAST.y_mc,i_mc );
    s=[0,cumsum(hypot(diff(COAST.x),diff(COAST.y)))];
    
    %% Force cyclic sections to be exactly cyclic
    cyclic=get_cyclic(COAST.x,COAST.y,COAST.ds0);
    if cyclic
       x1=.5*(COAST.x(1)+COAST.x(end));
       y1=.5*(COAST.y(1)+COAST.y(end));
       COAST.x(1)=x1;
       COAST.y(1)=y1;
       COAST.x(end)=x1;
       COAST.y(end)=y1;
       [COAST.x_mc,COAST.y_mc]=insert_section(COAST.x,COAST.y,COAST.x_mc,COAST.y_mc,i_mc);
    end
    [xx,yy,indi,indj]=get_intersections(COAST.x,COAST.y);
    xx0=unique(xx);
    
    merged=0;
    if length(xx0)==2
        indi=sort([indi(1:2),indj(1:2)]);
        % xnew=[[COAST.x(1:floor(indi(1))),xx(1),COAST.x(ceil(indi(4)):length(COAST.x))],nan,[xx(2),COAST.x(ceil(indi(2)):floor(indi(3))),xx(2)]];
        % ynew=[[COAST.y(1:floor(indi(1))),yy(1),COAST.y(ceil(indi(4)):length(COAST.y))],nan,[yy(2),COAST.y(ceil(indi(2)):floor(indi(3))),yy(2)]];
        xnew=[[COAST.x(1:floor(indi(1))),COAST.x(ceil(indi(4)):length(COAST.x))],nan,[COAST.x(ceil(indi(2)):floor(indi(3))),COAST.x(ceil(indi(2)))]];
        ynew=[[COAST.y(1:floor(indi(1))),COAST.y(ceil(indi(4)):length(COAST.y))],nan,[COAST.y(ceil(indi(2)):floor(indi(3))),COAST.y(ceil(indi(2)))]];
        merged=1;
        % insert new section in x_mc and y_mc
        [COAST.x_mc,COAST.y_mc]=insert_section(xnew,ynew,COAST.x_mc,COAST.y_mc,i_mc);
    elseif length(xx0)==1 
        indi=sort([indi(1),indj(1)]);
        xnew=[[COAST.x(1:floor(indi(1))),COAST.x(ceil(indi(2)):length(COAST.x))],nan,[COAST.x(ceil(indi(1)):floor(indi(2))),COAST.x(ceil(indi(1)))]];
        ynew=[[COAST.y(1:floor(indi(1))),COAST.y(ceil(indi(2)):length(COAST.y))],nan,[COAST.y(ceil(indi(1)):floor(indi(2))),COAST.y(ceil(indi(1)))]];
        merged=1;
        % insert new section in x_mc and y_mc
        [COAST.x_mc,COAST.y_mc]=insert_section(xnew,ynew,COAST.x_mc,COAST.y_mc,i_mc);
    end
end

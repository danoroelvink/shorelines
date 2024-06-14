function [COAST]=get_transportpoints(COAST,i_mc_range)
% function [COAST]=get_transportpoints(COAST,i_mc_range)
% 
% INPUT:
%   COAST
%        .x_mc          x-coordinate of coastal segment [m]
%        .y_mc          y-coordinate of coastal segment [m]
%        .ds0           default grid cell size
%        .BNDgroyne     index showing whether coastline element has a groyne at the beginning or end (.BNDgroyne=[0,1] means no groyne at the left and a groyne at the right side)
%   i_mc_range          index of the to be evaluated coastline segment/element, leave empty if all coastline elements should be evaluated.
%
% OUTPUT:
%   xq_mc               x-coordinate of transport points [m]
%   yq_mc               y-coordinate of transport points [m]
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2021 IHE Delft & Deltares
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

    % the function can be used for all coastline segments, or just this one
    if nargin==1
        i_mc_range=[1:COAST.n_mc];
        COAST.xq_mc=[];
        COAST.yq_mc=[];
    end

    % generate the transport grid for each of the coastline segments
    for i_mc=i_mc_range
        [x,y,n_mc]=get_one_polygon(COAST.x_mc,COAST.y_mc,i_mc);
        
        if length(x)>1
            % compute locations of transport points
            % and add extra boundary point
            xq=(x(1:end-1)+x(2:end))/2;
            yq=(y(1:end-1)+y(2:end))/2;
            cyclic=get_cyclic(x,y,COAST.ds0);
            if cyclic
                xq=[xq(end),xq,xq(1)];
                yq=[yq(end),yq,yq(1)];
            else
                if COAST.BNDgroyne(i_mc,1)==0       % extend grid when it is an open coast without a structure
                    xq=[1.5*x(1)-0.5*x(2),xq];
                    yq=[1.5*y(1)-0.5*y(2),yq];
                else                                % use first coastline point in case a groyne is present
                    xq=[x(1),xq];
                    yq=[y(1),yq];
                end
                if COAST.BNDgroyne(i_mc,2)==0       % extend grid when it is an open coast without a structure
                    xq=[xq,1.5*x(end)-0.5*x(end-1)];
                    yq=[yq,1.5*y(end)-0.5*y(end-1)];
                else                                % use last coastline point in case a groyne is present
                    xq=[xq,x(end)];
                    yq=[yq,y(end)];
                end
            end
        else
            fprintf('error : Too short coastline element with length 0\n');
            xq=x;
            yq=y;
        end
        
        %% reconstruct the xq_mc and yq_mc
        if isempty(COAST.xq_mc)
            COAST.xq_mc=xq;
            COAST.yq_mc=yq;
        elseif length(i_mc_range)>1 
            COAST.xq_mc=[COAST.xq_mc,nan,xq];
            COAST.yq_mc=[COAST.yq_mc,nan,yq];
        else
            [COAST.xq_mc,COAST.yq_mc]=insert_section(xq,yq,COAST.xq_mc,COAST.yq_mc,i_mc);
        end
    end
end



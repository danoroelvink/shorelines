function [ COAST,merged ] = merge_coastlines( COAST, i_mc)
% function [ COAST,merged ] = merge_coastlines( COAST, i_mc)
% 
% UNTITLED Summary of this function goes here
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
    [xx,yy]=get_intersections(COAST.x,COAST.y);
    if yesplot
        figure(3);
        plot(COAST.x,COAST.y,'.-b',xx,yy,'ok');
        hold on
        num=[1:length(COAST.x)];
        for i=1:length(COAST.x)
            text(COAST.x(i),COAST.y(i),num2str(num(i)));
        end
    end
    iX=0;
    ind=[];
    inp=[];
    for i=1:length(s)-1
        if ~isnan(s(i))&&~isnan(s(i+1))
            for ip=1:length(xx)
                err=abs(s(i+1)-s(i)-hypot(xx(ip)-COAST.x(i)  ,yy(ip)-COAST.y(i)) ...
                    -hypot(xx(ip)-COAST.x(i+1),yy(ip)-COAST.y(i+1)));
                if err<eps
                    iX=iX+1;
                    ind(iX)=i;
                    inp(iX)=ip;
                end
            end
        end
    end
    
    if isempty(ind) || length(ind)<2
        xnew=COAST.x;
        ynew=COAST.y;
        merged=0;
    elseif length(ind)<4
        %xnew=[COAST.x(1:ind(1)),COAST.x(ind(2)+1:end),nan,COAST.x(ind(1)+1:ind(2)),COAST.x(ind(1)+1)];
        %ynew=[COAST.y(1:ind(1)),COAST.y(ind(2)+1:end),nan,COAST.y(ind(1)+1:ind(2)),COAST.y(ind(1)+1)];
        xnew=[COAST.x(1:ind(1)),xx(inp(1)),COAST.x(ind(2)+1:end),nan,COAST.x(ind(1)+1:ind(2)),xx(inp(1)),COAST.x(ind(1)+1)];
        ynew=[COAST.y(1:ind(1)),yy(inp(1)),COAST.y(ind(2)+1:end),nan,COAST.y(ind(1)+1:ind(2)),yy(inp(1)),COAST.y(ind(1)+1)];
        merged=0;
    else
        xnew=[COAST.x(1:ind(1)),COAST.x(ind(4)+1:end),nan,COAST.x(ind(2)+1:ind(3)),COAST.x(ind(2)+1)];
        ynew=[COAST.y(1:ind(1)),COAST.y(ind(4)+1:end),nan,COAST.y(ind(2)+1:ind(3)),COAST.y(ind(2)+1)];
        % xnew=[COAST.x(1:ind(1)),xx(inp(1)),COAST.x(ind(4)+1:end),nan,COAST.x(ind(2)+1:ind(3)),xx(inp(2)),COAST.x(ind(2)+1)];
        % ynew=[COAST.y(1:ind(1)),yy(inp(1)),COAST.y(ind(4)+1:end),nan,COAST.y(ind(2)+1:ind(3)),yy(inp(2)),COAST.y(ind(2)+1)];
        merged=1;
    end
    if yesplot
        plot(xnew,ynew,'k','linewidth',2)
        hold off
    end

    % insert new section in x_mc and y_mc
    [COAST.x_mc,COAST.y_mc]=insert_section(xnew,ynew,COAST.x_mc,COAST.y_mc,i_mc);
    
    %% make transport points xq_mc and yq_mc
    [COAST]=get_transportpoints(COAST,i_mc);
end

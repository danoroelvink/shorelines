function [ shadowS,elementnr ] = find_shadows_mc(xq,yq,x_mc,y_mc,PHI,PHItdp,PHIbr,distw)
%
% Computes the indices of locations in the shadow zone 
% based on coastline shape, location of structures and wave incidence angle.
%
% INPUT:
%         x         : x-coordinate of coastline (only current section)
%         y         : y-coordinate of coastline (only current section)
%         x_mc      : x-coordinate of coastline (all sections)
%         y_mc      : y-coordinate of coastline (all sections)
%         PHI       : offshore wave incidence angle ([1] or [1xN] in degrees North)
%         PHItdp    : wave incidence angle at depth of closure ([1] or [1xN] in degrees North, using offshore wave if CERC is used)
%         PHIbr     : wave incidence angle at breaking point ([1] or [1xN] in degrees North, using offshore wave if CERC is used)
%         hard      : switch for hard structures
%
% OUTPUT:
%         xS        : x-coordinate of QS-points
%         yS        : y-coordinate of QS-points
%         shadowS   : Index of cells which are in the shadow zone
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

    % set QS-grid length and offshore wave angle
    nq=length(xq);   
    if length(PHI)==1
        PHI=repmat(PHI,[1,nq]);
    end
    if nargin<6
        PHItdp=PHI;
    end
    if nargin<7
        PHIbr=PHItdp;
    end
    if nargin<8
        distw=[];
    end
    
    % determine the length of the offshore, breaker and surfzone part of the wave rays
    lenSURF=repmat(250,[1,nq]);
    lenBREAKER=repmat(50,[1,nq]);
    len=repmat(5*hypot(max(xq)-min(xq),max(yq)-min(yq)),[1,nq]);
    if ~isempty(distw)
        lenBREAKER=min(lenBREAKER,distw*0.2);  % breaker zone cannot be more than 20% of distance to nearest wave output station 'distw'
        lenSURF=min(lenSURF,distw*0.8);        % surfzone cannot be more than 80% of distance to nearest wave output station 'distw'
        len=max(distw-lenBREAKER-lenSURF,1);   % offshore distance of ray is remainder that is not the surfzone or breaker zone, with minimum width of 1
    end
    
    % check wave rays to detect shadowing
    if nq==0
        shadowS=[];
    else
        elementnr={};
        shadowS=false(1,nq);
        elementnr=cell(nq,1);
        for i=1:nq
            % using a line width a small bend at the point of breaking (using PHIbr) and nearshore depth of closure (using PHItdp)
            xw=[xq(i)+1*sind(PHIbr(i)),...
                xq(i)+lenBREAKER(i)*sind(PHIbr(i)),...
                xq(i)+lenBREAKER(i)*sind(PHIbr(i))+lenSURF(i)*sind(PHItdp(i)),...
                xq(i)+lenBREAKER(i)*sind(PHIbr(i))+lenSURF(i)*sind(PHItdp(i))+len(i)*sind(PHI(i))];
            yw=[yq(i)+1*cosd(PHIbr(i)),...
                yq(i)+lenBREAKER(i)*cosd(PHIbr(i)),...
                yq(i)+lenBREAKER(i)*cosd(PHIbr(i))+lenSURF(i)*cosd(PHItdp(i)),...
                yq(i)+lenBREAKER(i)*cosd(PHIbr(i))+lenSURF(i)*cosd(PHItdp(i))+len(i)*cosd(PHI(i))];
            
            % identify if the waves are blocked by any other section
            [xx1,yy1,indc,inds]=get_intersections(x_mc,y_mc,xw,yw);
            shadowS(i)=~isempty(xx1);
            
            % identify which coastal section (or coastal structure) is blocking the waves
            if ~isempty(indc)
                section0=[];
                for cc=1:length(indc)
                jjcross=sum(isnan(x_mc(1:floor(indc(cc)))))+1;         
                section0(cc,1)=jjcross(1);
                end
                elementnr{i}=unique(section0);
            end
            
            if 0
                figure(10);clf
                % plot(xq,yq,x_mc,y_mc,xw,yw,'.-',xx1,yy1,'ro','linewidth',2);
                plot(xq,yq,'b+');hold on;
                plot(x_mc,y_mc,'r.-');
                plot(xw,yw,'g.-');
                %plot(P1(1,:),P1(2,:),'ro','linewidth',1);
                %xlim([min(xq),max(xq)]);
                %ylim([min(yq),max(yq)]);
                axis equal
                drawnow
                %%pause
            end
        end
    end
end

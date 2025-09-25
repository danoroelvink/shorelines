function [TRANSP,WAVE]=transport_shadow_treat(COAST,STRUC,WAVE,TRANSP)
%function [TRANSP,WAVE]=transport_shadow_treat(COAST,STRUC,WAVE,TRANSP)
%
% Computes the indices of grid cells with a shadowing of the waves 
% due to the coastline itself or hard structures.
%
% INPUT :  
%   COAST
%         .x           : x-coordinate of coastline [m] (only current section)
%         .y           : y-coordinate of coastline [m] (only current section)
%         .xq          : x-coordinate at transport points [m] (only current section)
%         .yq          : y-coordinate at transport points [m] (only current section)
%         .x_mc        : x-coordinate of coastline [m] (all sections)
%         .y_mc        : y-coordinate of coastline [m] (all sections)
%         .wberm       : berm/beach width along the coast [m]
%         .PHIcxy      : angle of the shorenormal of the coastline [°]
%   STRUC
%         .diffraction : swicth for using diffraction (0/1)
%         .xhard       : x-coordinate of hard structures [m]
%         .yhard       : y-coordinate of hard structures [m]
%         .xrevet      : x-coordinates of revetments [m]
%         .yrevet      : y-coordinates of revetments [m]
%   WAVE
%         .PHItdp      : wave incidence angle at depth-of-closure [°N]
%         .dPHIo       : relative angle of offshore waves with respect to the coastline [°]
%         .dPHItdp     : relative angle at depth-of-closure with respect to the coastline [°]
%         .dPHIbr      : relative angle at point of breaking with respect to the coastline [°]
%         .dPHIcrit    : critical wave angle for stability at depth-of-closure [°]
%         .dist        : closests distance of each of the wave input stations to the coastline [m]
%         .diffraction : switch for using diffraction (0/1)
%         .diff        : index with transport points affected by diffraction [-]
%   TRANSP
%         .QS          : transport rates [m3/yr]
%         .trform      : transport formulation used (e.g. 'CERC')
%         .xsedlim     : x-coordinates of regions with limited sediment availability / hard layers [m]
%         .ysedlim     : y-coordinates of regions with limited sediment availability / hard layers [m]
%
% OUTPUT :
%   TRANSP
%         .xS          : x-coordinate of QS-points [m]
%         .yS          : y-coordinate of QS-points [m]
%         .shadowS     : index of cells which are in the TRANSP.shadow zone of other sections of the coast (QS-points)
%         .shadowS_h   : index of cells which are in the TRANSP.shadow zone of hard structures (QS-points)
%         .shadow      : index of cells which are in the TRANSP.shadow zone of other sections of the coast (xy-point)
%         .shadowc     : index of cells which are in the TRANSP.shadow zone of hard structures (xy-point)
%   WAVE
%         .PHIbr       : angle of waves at point of breaking [°N] (only for CERC set to same value as .PHItdp)
%         .dPHItdp_sh  : relative incoming wave angle at depth-of-closure [°]
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

    nw=size(WAVE.HStdp,1); % number of wave conditions taken along at this timestep (can be more than 1 in case of simultaneous wave conditions)
    nq=size(WAVE.HStdp,2); % number of transport points
    nc=length(COAST.x);    % number of coastline points

    QS0=TRANSP.QS;
    WAVE.dPHItdp_sh=WAVE.dPHItdp;
    shadowS_90deg=abs(WAVE.dPHItdp)>2*WAVE.dPHIcrit; 
    
    rev.x=[];
    rev.y=[];
    if ~isempty(STRUC.xrevet) 
        rev.x=[rev.x,nan,STRUC.xrevet];
        rev.y=[rev.y,nan,STRUC.yrevet];
    end
    if ~isempty(TRANSP.xsedlim) 
        rev.x=[rev.x,nan,TRANSP.xsedlim];
        rev.y=[rev.y,nan,TRANSP.ysedlim];
    end
    struc.x=[];
    struc.y=[];
    if ~isempty(STRUC.xhard) 
        struc.x=[nan,STRUC.xhard,rev.x];
        struc.y=[nan,STRUC.yhard,rev.y];
    end
    if ~isempty(struc.x)
        struc.x=struc.x(2:end);
        struc.y=struc.y(2:end);
    end
    if ~isempty(rev.x)
        rev.x=rev.x(2:end);
        rev.y=rev.y(2:end);
    end

    if strcmpi(TRANSP.trform,'CERC')
        WAVE.PHIbr=WAVE.PHItdp;
    end
    
    % also account for local refraction in the nearshore
    [ TRANSP.shadowS ] = find_shadows_mc(COAST.xq,COAST.yq,COAST.x_mc,COAST.y_mc,WAVE.PHIo,WAVE.PHItdp,WAVE.PHIbr,WAVE.dist);
    TRANSP.shadowS = TRANSP.shadowS | shadowS_90deg;
    
    if ~isempty(STRUC.xhard) && ~isempty(COAST.x) && STRUC.diffraction==0
        %% shadowing in case without diffraction
        [ TRANSP.shadowS_h ] = find_shadows_mc(COAST.xq,COAST.yq,STRUC.xhard,STRUC.yhard,WAVE.PHIo,WAVE.PHItdp,WAVE.PHIbr,WAVE.dist);
        TRANSP.shadow=max([TRANSP.shadowS(:,1:end-1);TRANSP.shadowS_h(:,1:end-1);TRANSP.shadowS(:,2:end);TRANSP.shadowS_h(:,2:end)],[],1);
        TRANSP.shadowS_rev=zeros(nw,nq);
        if ~STRUC.diffraction
            TRANSP.shadow=max(TRANSP.shadowS(:,1:end-1),TRANSP.shadowS(:,2:end));
        end
        WAVE.dPHItdp_sh(TRANSP.shadowS)=0;
        WAVE.dPHItdp_sh(TRANSP.shadowS_h)=0;
    else
        %% shadowing in case with diffraction
        TRANSP.shadowS_h=false(nw,nq);
        WAVE.dPHItdp_sh(TRANSP.shadowS)=0;
        TRANSP.shadow=max(TRANSP.shadowS(:,1:end-1),TRANSP.shadowS(:,2:end));
    end

    %% add revetment shadows
    TRANSP.shadowS_rev=zeros(nw,nq);
    if ~isempty(rev.x) 
        xrev=rev.x;
        yrev=rev.y;
        [PHItdp]=get_smoothdata(WAVE.PHItdp,'angle',3);
        [ TRANSP.shadowS_revq ] = find_shadows_mc(COAST.xq,COAST.yq,xrev,yrev,WAVE.PHIo,PHItdp,PHItdp,mean(COAST.ds0*5));
        idshadowrevneg=find(TRANSP.shadowS_rev & TRANSP.QS<0);
        idshadowrevpos=find(TRANSP.shadowS_rev & TRANSP.QS>=0);
        [PHItdp]=get_smoothdata(WAVE.PHItdp,'anglemean',3);
        [PHIo]=get_smoothdata(WAVE.PHIo,'anglemean',3);
        [ TRANSP.shadowS_rev ] = find_shadows_mc(COAST.x,COAST.y,xrev,yrev,PHIo,PHItdp,PHItdp,mean(COAST.ds0*5));
        % idshadowrevneg=find(TRANSP.shadowS_rev & TRANSP.QS<0);
        % idshadowrevpos=find(TRANSP.shadowS_rev & TRANSP.QS>=0);
        % TRANSP.QS([TRANSP.shadowS_rev,false(nw,1)])=0;
        % TRANSP.QS([false(nw,1),TRANSP.shadowS_rev])=0;
        % TRANSP.QS(idshadowrevneg)=max(TRANSP.QS(min(idshadowrevneg+1,nq)),TRANSP.QS(idshadowrevneg));
        % TRANSP.QS(idshadowrevpos)=min(TRANSP.QS(max(idshadowrevpos-1,1)),TRANSP.QS(idshadowrevpos));
    end

    %% shadow zone for dune flux 
    if ~isempty(struc.x) && ~isempty(COAST.x) && ~isempty(COAST.wberm)
        xdune = COAST.x - COAST.wberm.*sind(COAST.PHIcxy); 
        ydune = COAST.y - COAST.wberm.*cosd(COAST.PHIcxy); 
        PHIoAvg = get_smoothdata(WAVE.PHIo,'anglemean');
        PHItdpAvg = get_smoothdata(WAVE.PHItdp,'anglemean');
        PHIbrAvg = get_smoothdata(WAVE.PHIbr,'anglemean');
        distAvg = (WAVE.dist(:,1:end-1)+WAVE.dist(:,2:end))/2;
        [ TRANSP.shadowS_hD ] = find_shadows_mc(xdune,ydune,struc.x,struc.y,PHIoAvg,PHItdpAvg,PHIbrAvg,distAvg);
    else
        TRANSP.shadowS_hD=false(nw,nc);
    end
    %%

    % Remove erroneous local points with just 1 shadowed cell
    dsh=diff(TRANSP.shadowS,1,2);
    idsh=find((dsh(:,1:end-1)-dsh(:,2:end))==2);
    if ~isempty(idsh)
        TRANSP.shadowS(idsh+1)=0;
    end
    
    % Exclude points with diffraction from the shadowing
    TRANSP.shadowS = TRANSP.shadowS & ~WAVE.diff; 
    TRANSP.shadowS_h = TRANSP.shadowS_h & ~WAVE.diff; 
    
    % SET QS=0 WHEN SHADOWING IS PRESENT
    usesmoothshadows=0;
    if usesmoothshadows==0
        TRANSP.QS(TRANSP.shadowS)=0;
        if length(TRANSP.shadowS_h)>0 && STRUC.diffraction==0
            TRANSP.QS(TRANSP.shadowS_h)=0;
        end
        
        % make sure sediment is not transported away from shadowed coastline points (by not allowing transport away from these points)
        if STRUC.diffraction==0
            for kk=1:nw
                idpos=find(TRANSP.QS(kk,:)>0 & [0,TRANSP.shadow(kk,:)]);
                TRANSP.QS(kk,idpos)=max(TRANSP.QS(kk,idpos-1),0);
                idneg=find(TRANSP.QS(kk,:)<0 & [TRANSP.shadow(kk,:),0]);
                TRANSP.QS(kk,idneg)=min(TRANSP.QS(kk,idneg+1),0);
            end
        end
        
    elseif usesmoothshadows==1
        pwr=0.5;                  % a higher factor will result in a more abrupt transition from shadowing to unshadowed (i.e. in regards to increase or decrease of QS)
        gaussianlengthfactor=3;   % Use 3 to extend influence to 1st grid cell at both sides next to transition. Use 5 to extend it to 2 grid cells beyond transition.
        TRANSP.shadowSfac1=1-get_smoothdata(TRANSP.shadowS,'gaussian', gaussianlengthfactor-1).^pwr;
        TRANSP.shadowSfac2=fliplr(1-get_smoothdata(fliplr(TRANSP.shadowS),'gaussian', gaussianlengthfactor-1).^pwr);
        TRANSP.shadowSfac=TRANSP.shadowSfac1;
        TRANSP.shadowSfac(TRANSP.QS<0)=TRANSP.shadowSfac2(TRANSP.QS<0);
        TRANSP.QS=TRANSP.shadowSfac.*TRANSP.QS;
        if length(TRANSP.shadowS_h)>0 && STRUC.diffraction==0 % to be tested with structure case & in combination with diffraction
            %TRANSP.QS(TRANSP.shadowS_h)=0;    % <- this was commented out in previous version
            TRANSP.shadowS_hfac=1-get_smoothdata(TRANSP.shadowS_h,'gaussian', 3);
            TRANSP.QS=TRANSP.shadowS_hfac.*TRANSP.QS;  
        end
    end
    % if breakerlinewaves==1
    %     dPHIbr_sh=dPHIbr;
    %     dPHIbr_sh(TRANSP.shadowS)=0;
    %     dPHIbr_sh(TRANSP.shadowS_h)=0;
    % end
    TRANSP.debug.QS1=TRANSP.QS;
    TRANSP.shadow=TRANSP.shadowS(:,1:end-1)&TRANSP.shadowS(:,2:end);             % this uses & instead of the | at line 89
    TRANSP.shadow_h=TRANSP.shadowS_h(:,1:end-1)&TRANSP.shadowS_h(:,2:end);
    
end

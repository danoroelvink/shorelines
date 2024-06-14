function [TRANSP,WAVE]=transport_shadow_treat(COAST,STRUC,WAVE,TRANSP)
%function [TRANSP,WAVE]=transport_shadow_treat(COAST,STRUC,WAVE,TRANSP)
%
% Computes the indices of grid cells with a shadowing of the waves 
% due to the coastline itself or hard structures.
%
% INPUT : 
%   COAST
%         .x         : x-coordinate of coastline (only current section)
%         .y         : y-coordinate of coastline (only current section)
%         .x_mc      : x-coordinate of coastline (all sections)
%         .y_mc      : y-coordinate of coastline (all sections)
%   STRUC
%         .x_hard    : x-coordinate of hard structures
%         .y_hard    : y-coordinate of hard structures
%   WAVE
%         .PHI       : Wave incidence angle ([1] or [1xN] in Radians and degrees North)
%         .dPHI      : Relative angle of offshore waves with respect to the coastline 
%         .diffraction
%
% OUTPUT:
%   TRANSP
%         .xS        : x-coordinate of QS-points
%         .yS        : y-coordinate of QS-points
%         .shadowS   : Index of cells which are in the TRANSP.shadow zone of other sections of the coast (QS-points)
%         .shadowS_h : Index of cells which are in the TRANSP.shadow zone of hard structures (QS-points)
%         .shadow    : Index of cells which are in the TRANSP.shadow zone of other sections of the coast (xy-point)
%         .shadowc   : Index of cells which are in the TRANSP.shadow zone of hard structures (xy-point)
%    dPHIcor   : Relative angle of offshore waves with respect to the coastline corrected to 0 in TRANSP.shadow zones
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

    QS0=TRANSP.QS;
    WAVE.dPHItdp_sh=WAVE.dPHItdp;
    shadowS_90deg=abs(WAVE.dPHItdp)>2*WAVE.dPHIcrit; 
    
    struc.x=[];
    struc.y=[];
    if ~isempty(STRUC.x_hard) 
        struc.x=[struc.x,nan,STRUC.x_hard];
        struc.y=[struc.y,nan,STRUC.y_hard];
    end
    if ~isempty(STRUC.x_revet) 
        struc.x=[struc.x,nan,STRUC.x_revet];
        struc.y=[struc.y,nan,STRUC.y_revet];
    end
    if ~isempty(TRANSP.x_sedlim) 
        struc.x=[struc.x,nan,TRANSP.x_sedlim];
        struc.y=[struc.y,nan,TRANSP.y_sedlim];
    end
    if ~isempty(struc.x)
        struc.x=struc.x(2:end);
        struc.y=struc.y(2:end);
    end
    
    if strcmpi(TRANSP.trform,'CERC')
        WAVE.PHIbr=WAVE.PHItdp;
    end
    
    % also account for local refraction in the nearshore
    [ TRANSP.shadowS ] = find_shadows_mc(COAST.xq,COAST.yq,COAST.x_mc,COAST.y_mc,WAVE.PHIo,WAVE.PHItdp,WAVE.PHIbr,WAVE.dist);
    TRANSP.shadowS = TRANSP.shadowS | shadowS_90deg;
    
    if ~isempty(struc.x) && ~isempty(COAST.x) && STRUC.diffraction==0
        [ TRANSP.shadowS_h ] = find_shadows_mc(COAST.xq,COAST.yq,struc.x,struc.y,WAVE.PHIo,WAVE.PHItdp,WAVE.PHIbr,WAVE.dist);
        TRANSP.shadow=max([TRANSP.shadowS(1:end-1);TRANSP.shadowS_h(1:end-1);TRANSP.shadowS(2:end);TRANSP.shadowS_h(2:end)],[],1);
    
        if ~STRUC.diffraction
            TRANSP.shadow=max(TRANSP.shadowS(1:end-1),TRANSP.shadowS(2:end));
        end
        WAVE.dPHItdp_sh(TRANSP.shadowS)=0;
        WAVE.dPHItdp_sh(TRANSP.shadowS_h)=0;
    else
        TRANSP.shadowS_h=false(1,length(TRANSP.QS));
        WAVE.dPHItdp_sh(TRANSP.shadowS)=0;
        TRANSP.shadow=max(TRANSP.shadowS(1:end-1),TRANSP.shadowS(2:end));
    end

    %% shadow zone for dune flux 
    if ~isempty(STRUC.x_revet) && ~isempty(COAST.x) && ~isempty(COAST.Wberm)
        xdune = COAST.x - COAST.Wberm.*sind(COAST.PHIcxy); 
        ydune = COAST.y - COAST.Wberm.*cosd(COAST.PHIcxy); 
        PHIoAvg = get_smoothdata(WAVE.PHIo,'anglemean');
        PHItdpAvg = get_smoothdata(WAVE.PHItdp,'anglemean');
        PHIbrAvg = get_smoothdata(WAVE.PHIbr,'anglemean');
        distAvg = (WAVE.dist(1:end-1)+WAVE.dist(2:end))/2;
        [ TRANSP.shadowS_hD ] = find_shadows_mc(xdune,ydune,STRUC.x_revet,STRUC.y_revet,PHIoAvg,PHItdpAvg,PHIbrAvg,distAvg);
    else
        TRANSP.shadowS_hD=false(1,length(COAST.x));
    end
    %%

    % Remove erroneous local points with just 1 shadowed cell
    dsh=diff(TRANSP.shadowS);
    idsh=find((dsh(1:end-1)-dsh(2:end))==2);
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
            idpos=find(TRANSP.QS>0 & [0,TRANSP.shadow]);
            TRANSP.QS(idpos)=max(TRANSP.QS(idpos-1),0);
            idneg=find(TRANSP.QS<0 & [TRANSP.shadow,0]);
            TRANSP.QS(idneg)=min(TRANSP.QS(idneg+1),0);
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
    TRANSP.shadow=TRANSP.shadowS(1:end-1)&TRANSP.shadowS(2:end);             % this uses & instead of the | at line 89
    TRANSP.shadow_h=TRANSP.shadowS_h(1:end-1)&TRANSP.shadowS_h(2:end);
    
end

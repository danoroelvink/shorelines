function [WAVE,TRANSP]=wave_angles(COAST,WAVE,TIDE,TRANSP,STRUC)
% function [WAVE,TRANSP]=wave_angles(COAST,WAVE,TIDE,TRANSP,STRUC)
%
% Computation of critical wave angles and QSmax at dynamic boundary and point of breaking.
% S-Phi curve at point of breaking to obtain 'critical coastline orientation' with respect to offshore conditions ('dPHI')
%
% INPUT: 
%    COAST
%       .nq                : number of qs-points
%    WAVE
%       .HStdp             : incoming wave height at depth-of-closure [m]
%       .PHItdp            : incoming wave direction at depth-of-closure [°N]
%       .TP                : wave period [s]
%       .dnearshore        : depth of dynamic boundary
%       .c1                : coefficient describing the curvature of the QS-Phi curve (only used in case trform='RAY')
%       .c2                : coefficient describing the curvature of the QS-Phi curve (only used in case trform='RAY')
%       .QSoffset          : coefficient describing the curvature of the QS-Phi curve (only used in case trform='RAY')
%    TRANSP 
%       .trform            : transport formulation (either 'CERC', 'KAMP', 'MILH', 'CERC3', 'VR14')
%       .gamma             : depth-induced breaking coefficient [-]
%       .alpha             : calibration factor for point of braking [-]
%       .qscal             : calibration coefficient for the transport computations [-]
%       .b                 : transport factor of CERC formula (only CERC)
%       .d50               : median grain size [um]
%       .d90               : 90-th percentile grain size [um]
%       .ks                : roughness height [m] 
%       .tanbeta           : slope [1:n] (only Kamphuis and Milhomens)
%       .pswell            : relative part of the pswell (in percentage, specify a value between 0 and 100)
%       .rhos              : density of sediment [kg/m3]
%       .rhow              : density of water [kg/m3]
%       .porosity          : porosity [-]
%       .g                 : acceleration of gravity [m/s2]
%       .suppresshighangle : switch to suppress high-angle instabilities
% 
% OUTPUT:
%    WAVE
%       .dPHIcrit          : Critical orientation of the coastline with respect to waves at the depth-of-closure (Nx2) [°]
%       .dPHIbr            : Critical orientation of the coastline with respect to waves at the point of breaking (Nx2) [°]
%    TRANSP
%       .QSmax             : Maximum transport foir critical wave angle [m3/yr]
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

    if ~strcmpi(TRANSP.trform,'RAY')
        PHIlocal = WAVE.PHItdp;
        HSlocal  = WAVE.HStdp;
        DEPlocal = repmat(WAVE.dnearshore,[1,COAST.nq]);
        [~,WAVE.ctdp,~,WAVE.ntdp]=get_disper(DEPlocal,WAVE.TP);
        if length(WAVE.TP)==1
            WAVE.TP=repmat(WAVE.TP,[1,COAST.nq]);
        end
    else
        c1 = WAVE.c1;
        c2 = WAVE.c2;
        QSoffset = WAVE.QSoffset;
    end
   
    %% Define critical angles for CERC and CERC2
    % critical angle for CERC1 = 45°    (Fixed, no influence of dynamic boundary)
    % critical angle for CERC2 = 42.39° (Fixed, no influence of dynamic boundary)
    % critical angle for other transport formulations depends on refraction and shoaling on the shoreface
    % these are therefore determined iteratively

    WAVE.dPHIcrit=repmat(45,[1,COAST.nq]);
    WAVE.dPHIcritbr=repmat(45,[1,COAST.nq]);
    WAVE0=WAVE;
    TRANSP.QSmax=zeros(1,COAST.nq);
    
    if strcmpi(TRANSP.trform,'CERC')
        WAVE.dPHIcrit=repmat(45,[1,COAST.nq]);
        WAVE0.dPHItdp=WAVE.dPHIcrit;
        TRANSP=transport(TRANSP,WAVE0,TIDE,STRUC,'QSmax');
        
    elseif strcmpi(TRANSP.trform,'CERC2')
        WAVE.dPHIcrit=repmat(42.39,[1,COAST.nq]);
        WAVE0.dPHItdp=WAVE.dPHIcrit;
        TRANSP=transport(TRANSP,WAVE0,TIDE,STRUC,'QSmax');
        
    elseif strcmpi(TRANSP.trform,'RAY')
        WAVE.dPHIcrit=abs(1./(sqrt(2).*c2));
        WAVE0.dPHItdp=WAVE.dPHIcrit;
        TRANSP=transport(TRANSP,WAVE0,TIDE,STRUC,'QSmax');
               
    elseif ~isempty(WAVE.sphimax)
        if isnumeric(WAVE.sphimax)
            % use a single sphimax value as estimate for the dPHIcrit
            WAVE.dPHIcrit=repmat(WAVE.sphimax,[1,COAST.nq]);
            WAVE.dPHIcritbr=repmat(WAVE.sphibrmax,[1,COAST.nq]);
            WAVE0.dPHItdp=WAVE.dPHIcrit;
            WAVE0.dPHIbr=WAVE.dPHIcritbr;
            TRANSP=transport(TRANSP,WAVE0,TIDE,STRUC,'QSmax');
        else % if strcmpi(WAVE.sphimax,'estimate')
            % initially compute the average dPHIcrit, and store it as a single sphimax value
            [WAVE,TRANSP]=get_Sphimax(WAVE,TIDE,TRANSP,STRUC);
            WAVE.sphimax=median(WAVE.dPHIcrit)-std(WAVE.dPHIcrit);
            WAVE.sphibrmax=median(WAVE.dPHIcritbr)-std(WAVE.dPHIcritbr);
        end
    else
        % compute dPHIcrit, dPHIcritbr and QSmax
        [WAVE,TRANSP]=get_Sphimax(WAVE,TIDE,TRANSP,STRUC);
    end
    
    %% Limit high-angle instabilities by forcing maximum transport when the angle is larger than dPHIcrit
    if TRANSP.suppresshighangle==1 
        for kk=1:size(WAVE.dPHItdp,1)
            id=abs(WAVE.dPHItdp(kk,:))>WAVE.dPHIcrit;
            WAVE.dPHItdp(kk,:)=sign(WAVE.dPHItdp(kk,:)).*min([abs(WAVE.dPHItdp(kk,:));WAVE.dPHIcrit],[],1);        
            WAVE.dPHIbr(kk,id)=sign(WAVE.dPHIbr(kk,id)).*min([abs(WAVE.dPHIbr(kk,id));WAVE.dPHIcritbr(id)],[],1);        
        end
    end
    
    %% limit refraction
    refraclimit=180.0;
    WAVE.dPHItdp((WAVE.dPHItdp-WAVE.dPHIo)<-refraclimit)=-refraclimit;
    WAVE.dPHItdp((WAVE.dPHItdp-WAVE.dPHIo)>refraclimit)=refraclimit;
    refraclimit=180.0;
    WAVE.dPHIbr((WAVE.dPHIbr-WAVE.dPHIo)<-refraclimit)=-refraclimit;
    WAVE.dPHIbr((WAVE.dPHIbr-WAVE.dPHIo)>refraclimit)=refraclimit;
    
end

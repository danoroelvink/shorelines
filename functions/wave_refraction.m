function [WAVE]=wave_refraction(WAVE,COAST,TRANSP)
% function [WAVE]=wave_refraction(WAVE,COAST,TRANSP)
%
% The wave refraction from the offshore point (ddeep) to the nearshore point (dnearshore) 
% is computed. Three options are available for the wave conditions. 
%    Option 1 : directly use offshore information (e.g. for 'CERC')
%    Option 2 : use snellius for wave transformation with shoaling and refraction (default)
%    Option 3 : extract waves from a 2D wave model
%
% INPUT: 
%    WAVE
%         .HSo              : significant wave height [m]
%         .PHIo             : wave direction [°N]
%         .TP               : wave period [s]
%         .ddeep            : depth at point of wave climate [m] (option 2)
%         .dnearshore       : depth of closure [m] (option 2)
%         .waveinteraction : switch for using the 2D wave model (option 3)
%         .wavefile         : filename of the wave model matrix (option 3)
%         .surfwidth        : surfzone width [m]
%    COAST
%         .x                : x-coordinates of the coastline (option 3)
%         .y                : y-coordinates of the coastline (option 3)
%         .n                : number of grid cells of the coastline (option 3)
%         .nq               : number of grid cells of the qs-points (option 3)
%         .PHIf             : orientation of the lower shoreface [deg N] (option 2)
%    TRANSP
%         .trform           : transport formulation (option 2)
%
% OUTPUT:
%    WAVE
%         .HStdp            : wave height at nearshore location along the coast [m] 
%         .PHItdp           : wave direction at nearshore location along the coast [°N] 
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
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses>
%   --------------------------------------------------------------------

    eps=0; 
    if WAVE.waveinteraction 
        % obtain waves from nearshore waves in MAT-file
        [ xg,yg,Hg_all,phiwg_all,dirtab,Tptab,Hstab ] = get_wave_fields_from_mat( WAVE.wavefile );
        phiwg_all=mod(phiwg_all,360.);
        if length(Hstab)>1
           [Hg,PHIwg,iwtw,PHIwg_wr,Hg_wr]=get_interpolated_wavefield_Hs_dir_Tp(xg,yg,Hg_all,phiwg_all,WAVE.HSo,WAVE.PHIo,WAVE.TP,Hstab,dirtab,Tptab);
        else
           [Hg,PHIwg,iwtw,PHIwg_wr,Hg_wr]=get_interpolated_wavefield_dir_Tp(xg,yg,Hg_all,phiwg_all,WAVE.HSo,WAVE.PHIo,WAVE.TP,dirtab,Tptab);
        end
        WAVE.PHIo=PHIwg;
        x=COAST.x;
        y=COAST.y;
        n=COAST.n;
        nq=COAST.nq;
        surfwidth=WAVE.surfwidth;
        [H,PHIcd]=get_interpolated_waves_from_grid(x,y,n,nq,surfwidth,xg,yg,Hg,PHIwg);
        WAVE.PHItdp=PHIcd;
        WAVE.HStdp=H;
        
    elseif ~isempty(COAST.PHIf) && ~strcmpi(TRANSP.trform,'CERC')  
        % nearshore conditions using snellius (shoaling and refraction)
        % compute relative angle
        relANGLE=mod(WAVE.PHIo-COAST.PHIf+180,360)-180;
        
        % make sure cells with relative incoming angle larger than 90 degrees have maximum refraction
        IDoff=find(abs(relANGLE)>=89.999); 
        relANGLE(IDoff)=sign(relANGLE(IDoff)) .* 89.999; 
        
        % compute nearshore conditions using snellius (shoaling and refraction)
        [~,c_deep,~,n_deep]=get_disper(WAVE.ddeep,WAVE.TP);
        [~,c_tdp,~,n_tdp]=get_disper(WAVE.dnearshore,WAVE.TP);
        WAVE.PHItdp=mod(COAST.PHIf+real(asind(c_tdp./c_deep.*sind(relANGLE))),360);
        relANGLEtdp=mod(WAVE.PHItdp-COAST.PHIf+180,360)-180;
        WAVE.HStdp=WAVE.HSo.*max(real(sqrt(n_deep.*c_deep.*cosd(relANGLE)./(n_tdp .*c_tdp .*cosd(relANGLEtdp)))),eps);
        
    else 
        % use offshore conditions at the toe of the dynamic profile
        WAVE.PHItdp=WAVE.PHIo;
        WAVE.HStdp=WAVE.HSo;
    end
    
    if length(WAVE.PHItdp)<=1
        WAVE.PHItdp=repmat(WAVE.PHItdp,[1,COAST.nq]);
    end
    if length(WAVE.HStdp)<=1
        WAVE.HStdp=repmat(WAVE.HStdp,[1,COAST.nq]);
    end
end

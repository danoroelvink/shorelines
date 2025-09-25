function [DUNE]=prepare_dunes(S,COAST)
% function [DUNE]=prepare_dunes(S,COAST)
%
% Initializes the parameters used for the dunes, wind and water-levels. 
% It also prepares the DUNE data-structure.
% 
% OUTPUT:
%    DUNE
%        .ds0
%        .ldbdune        : dune input file (overrules the xdune, ydune, wberm, dfelev, dcelev and cs parameters)
%        .xdune          : x coordinates where wberm and dfelev are given
%        .ydune          : y coordinates where wberm and dfelev are given
%        .wberm          : Berm width (distance MSL to dune foot) (m)
%        .dfelev         : Dune foot elevation (m)
%        .dcelev         : Dune crest elevation (m)
%        .cs             : Erosion coefficient of the dunes during storms, which scales the rate of erosion
%        .A              : Overwash parameter A
%        .duneaw         : Coefficient (Bagnold, 1937)
%        .rhos           : Density of the sediment (kg/m3)
%        .rhoa           : Density of air (kg/m3)
%        .g              : Gravitational acceleration (m/s2)
%        .d50            : Median grain diameter (m)
%        .d50r           : Median reference grain size (m)
%        .k              : Von Karman's coefficient
%        .kw             : Empirical coefficient (Sherman et al. 2013)
%        .porosity       : porosity
%        .segmaw         : Empirical factor used for scaling impact of the fetch length
%        .maxslope       : The maximum slope angle of the dunes (1:slope). The dunefoot height is lowered if the beach gets too steep (preserving the max slope).
%        .dtdune         : Timestep of the dune evolution
%      <cohesive/till>
%        .cstill         : Erosion coefficient of the dunes during storms for dunes with consolidated till layers, which scales the rate of erosion
%        .xtill          : Width of sandy dune in front of cohesive dune
%        .perctill       : Percentage of till (0 to 100) of the cohesive dune, with the sand percentage being 100-perctill
%      <wind>
%        .wndfile        : wind input file (overrules the uz, z, phiwnd0)
%        .uz             : wind velocity [m/s]
%        .z              : wind measurement vertical height [m]
%        .phiwnd0        : wind direction [°N]
%      <water-levels>
%        .watfile        : surge water-level file (overrules swl0)
%        .swl0           : static surge water-level [m w.r.t. MSL]
%        .runupform      : runup formulation applied
%        .runupfactor    : tuning factor for runup, considering inaccuracies in runup formulations
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

    fprintf('  Prepare dunes \n');
    DUNE=struct;
    DUNE.used=S.dune;
    DUNE.ldbdune=S.ldbdune;
    DUNE.wndfile=S.wndfile;
    DUNE.uz=S.uz;
    DUNE.z=S.z;
    DUNE.phiwnd0=S.phiwnd0;
    DUNE.watfile=S.watfile;
    DUNE.swl0=S.swl0;
    DUNE.cs=S.cs;
    DUNE.cstill=S.cstill;
    DUNE.xtill=S.xtill;
    DUNE.perctill=S.perctill;
    DUNE.A=S.aoverwash;
    DUNE.ds0=S.ds0;
    DUNE.duneaw=S.duneaw;
    DUNE.rhos=S.rhos;
    DUNE.rhoa=S.rhoa;
    DUNE.g=S.g;
    DUNE.d50=S.d50;
    DUNE.d50r=S.d50r;
    DUNE.k=S.k;
    DUNE.kw=S.kw;
    DUNE.porosity=S.porosity;
    DUNE.segmaw=S.segmaw;
    DUNE.xdune=S.xdune;               % x coords where wberm and dfelev are given
    DUNE.ydune=S.ydune;               % x coords where wberm and dfelev are given
    DUNE.xdune0=S.xdune;
    DUNE.ydune0=S.ydune;
    DUNE.wberm=S.wberm;               % Berm width (m)
    DUNE.wberm_mc0=[];                % Berm width at start simulation (m)
    DUNE.dfelev=S.dfelev;             % Dune foot elevation (m)
    DUNE.dcelev=S.dcelev;             % Dune crest elevation (m)
    DUNE.runupform=S.runupform;       % runup formulation applied
    DUNE.runupfactor=S.runupfactor;   % tuning factor for runup, considering inaccuracies in runup formulations
    DUNE.maxslope=S.maxslope;         % The maximum slope angle (1:slope). The dunefoot height is lowered if the beach gets too steep (preserving the max slope).
    DUNE.dt=S.dtdune;
    DUNE.csmodel=S.csmodel;           % file for the input of the CS-model (if non-empty, then the CS-model is used for the dunes isntead of the regular dune model)
    DUNE.temp=S.temp;
    
    %% Read dune properties at a number of longshore locations
    if DUNE.used 
        if ~isempty(findstr(lower(S.ldbdune),'.dun'))
            % LBDdune contains [x,y,  Bf0,  Bm0]
            dune=load(S.ldbdune);
            DUNE.xdune  = dune(:,1);
            DUNE.ydune  = dune(:,2);
            DUNE.wberm  = dune(:,3); % Berm width (distance MSL to dune foot) (m)
            DUNE.dfelev = dune(:,4); % Dune foot elevation (m)
            DUNE.dcelev = dune(:,5); % Dune crest elevation (m)
            if size(dune,2)>5
                DUNE.cs = dune(:,6);
            else
                DUNE.cs = repmat(DUNE.cs,size(DUNE.xdune));
            end
            if size(dune,2)>6
                DUNE.cstill = dune(:,7);
                DUNE.xtill = dune(:,8);
                DUNE.perctill = dune(:,9);
            else
                DUNE.cstill = repmat(DUNE.cstill,size(DUNE.xdune));
                DUNE.xtill = repmat(DUNE.xtill,size(DUNE.xdune));
                DUNE.perctill = repmat(DUNE.perctill,size(DUNE.xdune));
            end
            DUNE.wberm_mc0 = DUNE.wberm;
            DUNE.xdune0 = DUNE.xdune;
            DUNE.ydune0 = DUNE.ydune;
        end 

        %% write logfile
        % struct2log(DUNE,'DUNE','a');

    end
end

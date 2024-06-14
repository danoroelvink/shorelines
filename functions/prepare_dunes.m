function [DUNE]=prepare_dunes(S)
% function [DUNE]=prepare_dunes(S)
% 
% INPUT
%    S
%
% OUTPUT
%    DUNE
%       .xdune            x-coordinates 
%       .ydune            y-coordinates 
%       .Wberm            Berm width (distance MSL to dune foot) (m)
%       .Dfelev           Dune foot elevation (m)
%       .Dcelev           Dune crest elevation (m)
%
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

    fprintf('  Prepare dunes \n');
    DUNE=struct;
    DUNE.used = S.dune;
    DUNE.LDBdune = S.LDBdune;
    DUNE.WNDfile = S.WNDfile;
    DUNE.uz = S.uz;
    DUNE.z = S.z;
    DUNE.phiwnd0 = S.phiwnd0;
    DUNE.WATfile = S.WATfile;
    DUNE.SWL0 = S.SWL0;
    DUNE.Cs = S.Cs;
    DUNE.Cstill = S.Cstill;
    DUNE.xtill = S.xtill;
    if isempty(DUNE.xtill)
        DUNE.xtill=inf;
    end
    DUNE.perctill = S.perctill;
    DUNE.A = S.A_overwash;
    DUNE.ds0 = S.ds0;
    DUNE.duneAw = S.duneAw;
    DUNE.rhos = S.rhos;
    DUNE.rhoa = S.rhoa;
    DUNE.g = S.g;
    DUNE.d50 = S.d50;
    DUNE.d50r = S.d50r;
    DUNE.k = S.k;
    DUNE.Kw = S.Kw;
    DUNE.porosity = S.porosity;
    DUNE.segmaw = S.segmaw;
    DUNE.xdune=[];                    % x coords where Wberm and Dfelev are given
    DUNE.ydune=[];                    % x coords where Wberm and Dfelev are given
    DUNE.Wberm=[];                    % Berm width (m)
    DUNE.Wberm_mc0=[];                % Berm width at start simulation (m)
    DUNE.Dfelev=[];                   % Dune foot elevation (m)
    DUNE.Dcelev=[];                   % Dune crest elevation (m)
    DUNE.runupform=S.runupform;       % runup formulation applied
    DUNE.runupfactor=S.runupfactor;   % tuning factor for runup, considering inaccuracies in runup formulations
    DUNE.maxslope=S.maxslope;         % The maximum slope angle (1:slope). The dunefoot height is lowered if the beach gets too steep (preserving the max slope).
    DUNE.dt=S.dtdune;
    
    %% Read dune properties at a number of longshore locations
    if DUNE.used
       if ~isempty(findstr(lower(S.LDBdune),'.dun'))
          % LBDdune contains [x,y,  Bf0,  Bm0]
          dune=load(S.LDBdune);
          DUNE.xdune  = dune(:,1);
          DUNE.ydune  = dune(:,2);
          DUNE.Wberm  = dune(:,3); % Berm width (distance MSL to dune foot) (m)
          DUNE.Dfelev = dune(:,4); % Dune foot elevation (m)
          DUNE.Dcelev = dune(:,5); % Dune crest elevation (m)
          if size(dune,2)>5
              DUNE.Cs = dune(:,6);
          else
              DUNE.Cs = repmat(DUNE.Cs,size(DUNE.xdune));
          end
          if size(dune,2)>6
              DUNE.Cstill = dune(:,7);
              DUNE.xtill = dune(:,8);
              DUNE.perctill = dune(:,9);
          else
              DUNE.Cstill = repmat(DUNE.Cstill,size(DUNE.xdune));
              DUNE.xtill = repmat(DUNE.xtill,size(DUNE.xdune));
              DUNE.perctill = repmat(DUNE.perctill,size(DUNE.xdune));
          end
          DUNE.Wberm_mc0 = DUNE.Wberm;
          DUNE.xdune0 = DUNE.xdune;
          DUNE.ydune0 = DUNE.ydune;
       end 

       %% write logfile
       % struct2log(DUNE,'DUNE','a');

    end

end

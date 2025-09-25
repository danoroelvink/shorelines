function [COAST]=interpolate_props(COAST,DUNE,MUD,TIME)
% function [COAST]=interpolate_props(COAST,DUNE,MUD,TIME)
% 
% This function re-interpolates parameters along the coast at 
% every timestep, which is needed when the grid changes. 
% 
% INPUT: 
%      COAST
%         .x          : x-coordinate of transport points (m)
%         .y          : y-coordinate of transport points (m)
%      DUNE
%         .xdune      : x coordinates where dune props are specified (m)
%         .ydune      : y coordinates where dune props are specified
%         .wberm      : initial beach width at these locations (m)
%         .dfelev     : initial dune foot elevation (m)
%         .dcelev     : initial dune crest elevation (m)
%         .wberm_mc0  : initial beach width at start of simulation (m)
%      MUD
%         .xmgv       : x coordinates where mangrove props are specified (m)
%         .ymgv       : y coordinates where mangrove props are specified (m)
%         .Bf         : initial mud flat width (m)
%         .Bm         : initial mangrove width (m)
%
% OUTPUT:
%      COAST
%         .wberm      : beach width
%         .dfelev     : dune foot elevation
%         .dcelev     : dune crest elevation
%         .Bf         : mud flat width
%         .Bm         : mangrove width
%         .Bfm        : Mangrove width on mud flat
%         .cs         : erosion factor of sandy dunes
%         .cstill     : erosion factor of cohesive/till dunes 
%         .xtill      : width of the sand dune in front of the cohesive dune
%         .tillperc   : percentage of cohesive material in the dunes (sand percentage = 100-perctill)
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
   if ~MUD.used && ~DUNE.used
      COAST.wberm = [];
      return;
   end
   
   i_mc=COAST.i_mc;
   if TIME.it>=0 %
      %% Interpolation of dune and mangrove properties for all timesteps
      % Initially, values read from file are used; in subsequent steps
      % update_props updates DUNE and MUD
      x            =   COAST.x;
      y            =   COAST.y;
      method       =   'weighted_distance';
      
      if DUNE.used
          % Interpolate sandy dune properties alongshore
          fieldnm1={'wberm','dfelev','dcelev'};
          [var1i,idGRID,distw]=get_interpolation_on_grid(method,x,y,DUNE.xdune,DUNE.ydune,DUNE,fieldnm1,{});
          COAST.wberm  = var1i.wberm(:)';
          COAST.dfelev  = var1i.dfelev(:)';
          COAST.dcelev  = var1i.dcelev(:)';
          
          % Interpolate cohesive dune properties alongshore
          fieldnm1={'cs'};
          if ~isempty(DUNE.xtill)
          fieldnm1={'cs','cstill','xtill','perctill'};
          end
          [var1i,idGRID,distw]=get_interpolation_on_grid(method,x,y,DUNE.xdune0,DUNE.ydune0,DUNE,fieldnm1,{});
          COAST.cs = var1i.cs(:)';
          if ~isempty(DUNE.xtill)
              COAST.cstill = var1i.cstill(:)';
              COAST.perctill = var1i.perctill(:)';
              COAST.xtill  = var1i.xtill(:)';
          else
              COAST.cstill=COAST.cs;
              COAST.perctill = zeros(size(COAST.cs));
              COAST.xtill  = [];
          end
          
          % Make sure to either initialize the 'xtill' or re-interpolate. 
          if TIME.it~=0 && ~isempty(DUNE.xtill)
              x1_mc=COAST.x1_mc;
              y1_mc=COAST.y1_mc;
              if COAST.i_mc==1
                  COAST.xtill_mc1=COAST.xtill_mc(:)'; % cross-shore position of till in dune [m w.r.t. dune front]
              end
              fieldnm1={'xtill_mc1'};         
              COAST.xtill_mc1=interpNANs(COAST.xtill_mc1);
              [var1i]=get_interpolation_on_grid(method,x,y,x1_mc,y1_mc,COAST,fieldnm1,{});
              COAST.xtill=var1i.xtill_mc1(:)';
          end
      end
      if MUD.used
          % Interpolate mangrove properties alongshore
          fieldnm1={'Bf','Bm','Bfm'};
          [var1i,~,~]=get_interpolation_on_grid(method,x,y,MUD.xmgv,MUD.ymgv,MUD,fieldnm1,{});
          for ff=1:length(fieldnm1)
              COAST.(fieldnm1{ff})  = var1i.(fieldnm1{ff})(:)';
          end
      end
   end
end

function [COAST]=interpolate_props(COAST,DUNE,MUD,TIME)
% function [COAST]=interpolate_props(COAST,DUNE,MUD,TIME)
%
% INPUT :
%      COAST
%         .x           : x-coordinate of transport points (m)
%         .y           : y-coordinate of transport points (m)
%      DUNE
%         .xdune       : x coordinates where dune props are specified (m)
%         .ydune       : y coordinates where dune props are specified
%         .Wberm       : initial beach width at these locations (m)
%         .Dfelev      : initial dune foot elevation (m)
%         .Dcelev      : initial dune crest elevation (m)
%         .Wberm_mc0   : initial beach width at start of simulation (m)
%      MUD
%         .xmgv        : x coordinates where mangrove props are specified (m)
%         .ymgv        : y coordinates where mangrove props are specified (m)
%         .Bf          : initial mud flat width (m)
%         .Bm          : initial mangrove width (m)
%
% OUTPUT :
%      COAST
%         .Wberm    : beach width
%         .Dfelev   : dune foot elevation
%         .Dcelev   : dune crest elevation
% 
%         .Bf       : mud flat width
%         .Bm       : mangrove width
%         .Bfm      : Mangrove width on mud flat
%
%         .Cs       : 
%         .Cstill   : 
%         .xtill    : 
%         .tillperc : 
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
   if ~MUD.used && ~DUNE.used
      COAST.Wberm = [];
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
          % fprintf('  Interpolate dune properties alongshore \n');
          %fldnms1={'Wberm','Dfelev','Dcelev','xhard','qshard','xtill','perctill'};
          fieldnm1={'Wberm','Dfelev','Dcelev'};
          [var1i,idGRID,distw]=get_interpolation_on_grid(method,x,y,DUNE.xdune,DUNE.ydune,DUNE,fieldnm1,{});
          COAST.Wberm  = var1i.Wberm(:)';
          COAST.Dfelev  = var1i.Dfelev(:)';
          COAST.Dcelev  = var1i.Dcelev(:)';
          fieldnm1={'Cs','Cstill','xtill','perctill'};
          [var1i,idGRID,distw]=get_interpolation_on_grid(method,x,y,DUNE.xdune0,DUNE.ydune0,DUNE,fieldnm1,{});
          COAST.Cs = var1i.Cs(:)';
          COAST.Cstill = var1i.Cstill(:)';
          COAST.perctill = var1i.perctill(:)';
          if TIME.it==0  %|| (~strcmpi(fieldnm1{ff},'xtill') && ~strcmpi(fieldnm1{ff},'xhard'))
              COAST.xtill  = var1i.xtill(:)';
              %COAST.xhard  = var1i.xhard(:)';
          else
              x1_mc=COAST.x1_mc;
              y1_mc=COAST.y1_mc;
              if COAST.i_mc==1
                  COAST.xtill_mc1=COAST.xtill_mc(:)'; % cross-shore position of till in dune [m w.r.t. dune front]
                  %COAST.xhard_mc1=COAST.xhard_mc(:)'; % cross-shore position of hard layer in coast [m w.r.t. coastline]
              end
              fieldnm1={'xtill_mc1'};             
              [var1i]=get_interpolation_on_grid(method,x,y,x1_mc,y1_mc,COAST,fieldnm1,{});
              COAST.xtill=var1i.xtill_mc1(:)';
          end
      end
      if MUD.used
          % fprintf('  Interpolate mangrove properties alongshore \n');
          fieldnm1={'Bf','Bm','Bfm'};
          [var1i,~,~]=get_interpolation_on_grid(method,x,y,MUD.xmgv,MUD.ymgv,MUD,fieldnm1,{});
          for ff=1:length(fieldnm1)
              COAST.(fieldnm1{ff})  = var1i.(fieldnm1{ff})(:)';
          end
      end
   end
end

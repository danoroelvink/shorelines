function [COAST,DUNE,MUD] = update_props(COAST,DUNE,MUD)
% function [COAST,DUNE,MUD] = update_props(COAST,DUNE,MUD)
% 
% Updates the dune and mud properties when the grid has changed. 
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
   if DUNE.used 
      ok = ~isnan(COAST.x_mc);
      DUNE.xdune = COAST.x_mc(ok); % - COAST.wberm_mc(ok).*sind(COAST.PHIcxy_mc(ok));
      DUNE.ydune = COAST.y_mc(ok); % - COAST.wberm_mc(ok).*cosd(COAST.PHIcxy_mc(ok));
      DUNE.wberm = COAST.wberm_mc(ok);
      DUNE.dfelev = COAST.dfelev_mc(ok);
      DUNE.dcelev = COAST.dcelev_mc(ok);
      COAST.xdune_mc = COAST.x_mc - COAST.wberm_mc.*sind(COAST.PHIcxy_mc); 
      COAST.ydune_mc = COAST.y_mc - COAST.wberm_mc.*cosd(COAST.PHIcxy_mc); 
      %[COAST]=get_disentangled(COAST,{'xdune_mc','ydune_mc'});
   end
   if MUD.used
      ok = ~isnan(COAST.x_mc);
      MUD.xmgv = COAST.x_mc(ok);
      MUD.ymgv = COAST.y_mc(ok);
      MUD.Bf = COAST.Bf_mc(ok);
      MUD.Bfm = COAST.Bfm_mc(ok);
      MUD.Bm = COAST.Bm_mc(ok);
   end

end
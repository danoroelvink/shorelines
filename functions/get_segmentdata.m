function [COAST,WAVE,TRANSP,i1,i2]=get_segmentdata(COAST,WAVE,TRANSP,DUNE,MUD,i_mc)
% function [COAST,WAVE,TRANSP,i1,i2]=get_segmentdata(COAST,WAVE,TRANSP,DUNE,MUD,i_mc)
% 
% The data of the active coastal element 'i_mc' is collected from the variables 
% that contain information from all coastal elements (the '_mc' elements)
% 
% INPUT: 
%    COAST
%      .x_mc           : x-coordinates for all coastal elements
%      .y_mc           : y-coordinates for all coastal elements
%      .s_mc           : length along the coast for all coastal elements
%      .ds_mc          : grid cell size for all cells of all coastal elements
%      .h0_mc          : active profile height for all coastal elements
%    WAVE
%      .HSo_mc         : wave height at offshore point [m]
%      .TP_mc          : wave period [s]
%      .PHItdp_mc      : incoming nearhsore wave angle for all coastal elements
%      .PHIbr_mc       : incoming wave angle at point of breaking for all coastal elements
%      .diff_mc        : indices of coastline points affected by wave diffraction
%    TRANSP
%      .trform         : transport formulation
%      .QS_mc          : transports for all coastal elements
%      .shadow_mc      : index of cells with shadowing of other parts of the coastlines (xy-points)
%      .shadowS_hD_mc  : index of cells with shadowing of hard structures (QS-points)
%    DUNE
%      .qs_mc          : dune erosion volume change that increases the beach width [m3/m/yr]
%      .qss_mc         : dune erosion volume change that increases the beach width, part that is sandy [m3/m/yr]
%      .ql_mc          : dune erosion volume change that does not increase the beach width [m3/m/yr]
%      .qw_mc          : wind transport from the beach to the dune [m3/m/yr]
%      .R_mc           : runup level for dune erosion [m]
%      .SWL_mc         : still water level for dune erosion [m]
%      .wberm_mc       : width of the berm/beach [m]
%      .dfelev_mc      : Dune foot elevation (m)
%      .dcelev_mc      : Dune crest elevation (m)
%    MUD
%      .dndt_mud_mc    : change in mud flat position in [m3/yr/m]
%      .dBfdt_mc       : change in mudflat width [m/yr]
%      .dBmdt_mc       : change in mangove width [m/yr]
%      .dBfmdt_mc      : change in colonizing mangrove width [m/yr] 
%      .Bf_mc          : mud flat width [m]
%      .Bm_mc          : mangrove width [m]
%      .Bfm_mc         : colonizing mangrove width [m]   
%    i_mc              : index of segment to be retrieved
%         
% OUTPUT:
%    COAST
%      .x              : x-coordinates of cells for considered coastal segment
%      .y              : y-coordinates of cells for considered coastal segment
%      .n              : number of coastline points of considered coastal segment
%      .nq             : number of qs-points of considered coastal segment
%      .s              : length along the coast for considered coastal segment
%      .ds             : grid cell size for all cells of considered segment
%      .h0             : active profile height of considered segment
%      .cyclic         : index whether considered coastal segment is cyclical
%      .i1             : start index of this element in x_mc
%      .i2             : end index of this element in x_mc
%      .i_mc           : number of active coastal element 
%      .n_mc           : number of coastal elements
%      .qs             : dune erosion volume change that increases the beach width [m3/m/yr]
%      .qss            : dune erosion volume change that increases the beach width, part that is sandy [m3/m/yr]
%      .ql             : dune erosion volume change that does not increase the beach width [m3/m/yr]
%      .qw             : wind transport from the beach to the dune [m3/m/yr]
%      .R              : runup level for dune erosion [m]
%      .SWL            : still water level for dune erosion [m]
%      .wberm          : width of the berm/beach [m]
%      .dfelev         : Dune foot elevation (m)
%      .dcelev         : Dune crest elevation (m)
%      .dndt_mud       : change in mud flat position in [m3/yr/m]
%      .dBfdt          : change in mudflat width [m/yr]
%      .dBmdt          : change in mangove width [m/yr]
%      .dBfmdt         : change in colonizing mangrove width [m/yr] 
%      .Bf             : mud flat width [m]
%      .Bm             : mangrove width [m]
%      .Bfm            : colonizing mangrove width [m]
%    WAVE
%      .HSo            : wave height at offshore point [m]
%      .TP             : wave period [s]
%      .PHI            : incoming wave angle (either nearshore or at point of breaking) for considered coastal segment
%      .diff           : indices of coastline points affected by wave diffraction
%    TRANSP
%      .QS             : transports for considered coastal segment
%      .shadow         : index of cells with shadowing of other parts of the coastlines (xy-points)
%      .shadowS_hD     : index of cells with shadowing of hard structures (QS-points)
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

    COAST.i_mc=i_mc; 
    [ COAST.x,COAST.y,COAST.n_mc,i1,i2 ] = get_one_polygon( COAST.x_mc,COAST.y_mc,i_mc );
    [ COAST.s ] = get_one_polygon(COAST.s_mc,i_mc);
    [ COAST.ds ] = get_one_polygon(COAST.ds_mc,i_mc);
    [ COAST.h0 ] = get_one_polygon(COAST.h0_mc,i_mc);
    [ COAST.PHIcxy ] = get_one_polygon(COAST.PHIcxy_mc,i_mc);
    [ COAST.PHIcxy0 ] = get_one_polygon(COAST.PHIcxy0_mc,i_mc);
    [ WAVE.HSo ] = get_one_polygon( WAVE.HSo_mc,i_mc);
    [ WAVE.TP ] = get_one_polygon( WAVE.TP_mc,i_mc);
    if DUNE.used
        [ COAST.qs ]       = get_one_polygon( COAST.qs_mc,i_mc);
        [ COAST.qss ]      = get_one_polygon( COAST.qss_mc,i_mc);
        [ COAST.ql ]       = get_one_polygon( COAST.ql_mc,i_mc);
        [ COAST.qw ]       = get_one_polygon( COAST.qw_mc,i_mc);
        [ COAST.R ]        = get_one_polygon( COAST.R_mc,i_mc);
        [ COAST.SWL ]      = get_one_polygon( COAST.SWL_mc,i_mc);
        [ COAST.wberm ]    = get_one_polygon( COAST.wberm_mc,i_mc);
        [ COAST.dfelev ]   = get_one_polygon( COAST.dfelev_mc,i_mc);
        [ COAST.dcelev ]   = get_one_polygon( COAST.dcelev_mc,i_mc);
    end
    if MUD.used
        [ COAST.dndt_mud ] = get_one_polygon( COAST.dndt_mudmc,i_mc);
        [ COAST.dBfdt ]    = get_one_polygon( COAST.dBfdt_mc,i_mc);
        [ COAST.dBmdt ]    = get_one_polygon( COAST.dBmdt_mc,i_mc);
        [ COAST.dBfmdt ]   = get_one_polygon( COAST.dBfmdt_mc,i_mc);
        [ COAST.Bf ]       = get_one_polygon( COAST.Bf_mc,i_mc);
        [ COAST.Bm ]       = get_one_polygon( COAST.Bm_mc,i_mc);
        [ COAST.Bfm ]      = get_one_polygon( COAST.Bfm_mc,i_mc);
    end
    if (strcmpi(TRANSP.trform,'CERC') || strcmpi(TRANSP.trform,'CERC2'))
        [ TRANSP.QS,WAVE.PHI ] = get_one_polygon(TRANSP.QS_mc,WAVE.PHItdp_mc,i_mc);
    else
        [ TRANSP.QS,WAVE.PHI ] = get_one_polygon(TRANSP.QS_mc,WAVE.PHIbr_mc,i_mc);
    end
    COAST.shadow_mc = TRANSP.shadow_mc;
    [ TRANSP.shadow ] = get_one_polygon(TRANSP.shadow_mc,i_mc);
    [ TRANSP.shadowS_hD ] = get_one_polygon(TRANSP.shadowS_hD_mc,i_mc);
    COAST.shadow=TRANSP.shadow;
    [ ~, WAVE.diff ] = get_one_polygon(WAVE.PHIbr_mc,WAVE.diff_mc,i_mc);
    COAST.cyclic = get_cyclic(COAST.x,COAST.y,COAST.ds0);
    COAST.n=length(COAST.x);
    COAST.nq=length(TRANSP.QS);
    COAST.i_mc=i_mc;
    COAST.i1=i1;
    COAST.i2=i2;
    
end

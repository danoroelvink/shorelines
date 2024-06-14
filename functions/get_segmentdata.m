function [COAST,WAVE,TRANSP,i1,i2]=get_segmentdata(COAST,WAVE,TRANSP,DUNE,MUD,i_mc)
% [COAST,WAVE,TRANSP,i1,i2]=get_segmentdata(COAST,WAVE,TRANSP,DUNE,MUD,i_mc)
%
% INPUT:
%    COAST
%      .x_mc               x-coordinates for all coastal segments
%      .y_mc               y-coordinates for all coastal segments
%      .h0_mc              Active profile height for all coastal segments
%    WAVE
%      .PHItdp_mc          incoming nearhsore wave angle for all coastal segments
%      .PHIbr_mc           incoming wave angle at point of breaking for all coastal segments
%    TRANSP
%      .QS_mc              transports for all coastal segments
%      .trform             transport formulation
%    i_mc                  index of segment to be retrieved
%         
% OUTPUT:
%    COAST
%      .x                  x-coordinates of cells for considered coastal segment
%      .y                  y-coordinates of cells for considered coastal segment
%      .n                  number of cells for considered coastal segment
%      .s                  length along the coast for considered coastal segment
%      .cyclic             index whether considered coastal segment is cyclical
%      .n_mc               number of coastal segments
%    WAVE
%      .PHI                incoming wave angle (either nearshore or at point of breaking) for considered coastal segment
%    TRANSP
%      .QS                 transports for considered coastal segment
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

    COAST.i_mc=i_mc; 
    [ COAST.x,COAST.y,COAST.n_mc,i1,i2 ] = get_one_polygon( COAST.x_mc,COAST.y_mc,i_mc );
    [ COAST.s ] = get_one_polygon(COAST.s_mc,i_mc);
    [ COAST.ds ] = get_one_polygon(COAST.ds_mc,i_mc);
    [ COAST.h0 ] = get_one_polygon(COAST.h0_mc,i_mc);
    [ WAVE.HSo ] = get_one_polygon( WAVE.HSo_mc,i_mc);
    [ WAVE.TP ] = get_one_polygon( WAVE.TP_mc,i_mc);
    if DUNE.used
        [ COAST.qs ]       = get_one_polygon( COAST.qs_mc,i_mc);
        [ COAST.ql ]       = get_one_polygon( COAST.ql_mc,i_mc);
        [ COAST.qw ]       = get_one_polygon( COAST.qw_mc,i_mc);
        [ COAST.R ]        = get_one_polygon( COAST.R_mc,i_mc);
        [ COAST.SWL ]      = get_one_polygon( COAST.SWL_mc,i_mc);
        [ COAST.Wberm ]    = get_one_polygon( COAST.Wberm_mc,i_mc);
        [ COAST.Dfelev ]   = get_one_polygon( COAST.Dfelev_mc,i_mc);
        [ COAST.Dcelev ]   = get_one_polygon( COAST.Dcelev_mc,i_mc);
    end
    if MUD.used
        [ COAST.dndt_mud ] = get_one_polygon( COAST.dndt_mud_mc,i_mc);
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

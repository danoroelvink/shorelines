function [SCALARS,VECTORS,idgrid]=get_gridprojection(COAST,WAVE,TRANSP,xp,yp,xw,yw,scalarfields,vectorfields,method)
% function [SCALARS,VECTORS,idgrid]=get_gridprojection(COAST,WAVE,TRANSP,xp,yp,xw,yw,scalarfields,vectorfields,method)
%
% The routine converts a variable from a coastline to another interpolated grid.
% 
% INPUT: 
%      COAST          : structure with coast info, and fields
%      WAVE           : structure with wave info, and fields
%      TRANSP         : structure with transport information, and fields
%      scalarfields   : parameters to be projected, of fields, e.g. {'HStdp'}
%      vectorfields   : parameters to be projected, of fields, e.g. {'PHItdp'}
%      method         : interpolation method 'weighted_distance' or 'alongshore_mapping'
% 
% OUTPUT:
%      scalars        : interpolated scalars
%      vectors        : interpolated vectors
%      idgrid         : interpolation grid indices
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

    % get scalar fields
    var1=struct;
    for ff=1:length(scalarfields)
        if isfield(WAVE,[scalarfields{ff},'_mc'])
            var1.(scalarfields{ff})=WAVE.([scalarfields{ff},'_mc']);
        elseif isfield(TRANSP,[scalarfields{ff},'_mc'])
            var1.(scalarfields{ff})=TRANSP.([scalarfields{ff},'_mc']);
        else
            var1.(scalarfields{ff})=COAST.([scalarfields{ff},'_mc']);
        end
    end
    %var1.h0=get_convertgrid(var1.h0,GROYNE,COAST);
    
    % get vector fields
    var2=struct;
    for ff=1:length(vectorfields)
        if isfield(WAVE,[vectorfields{ff},'_mc'])
            var2.(vectorfields{ff})=WAVE.([vectorfields{ff},'_mc']);
        elseif isfield(TRANSP,[vectorfields{ff},'_mc'])
            var2.(vectorfields{ff})=TRANSP.([vectorfields{ff},'_mc']);
        else
            var2.(vectorfields{ff})=COAST.([vectorfields{ff},'_mc']);
        end
    end
    
    % use weighted distance to interpolate info on the grid
    [SCALARS,VECTORS,idgrid]=get_interpolation_on_grid(method,xp,yp,xw,yw,var1,var2);
end
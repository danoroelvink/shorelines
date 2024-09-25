function [COAST] = get_active_profile(COAST)
% function [COAST] = get_active_profile(COAST)
% 
% This routine obtains the spatially varying or static active height of the cross-shore profile.
%
% INPUT:
%    COAST     : Structure with data of the ShorelineS model, of which is used
%      .x      : x-coordinate of coastal segment [m]
%      .y      : y-coordinate of coastal segment [m]
%      .n      : number of grid cells of coastal segment
%      .h0     : Active profile height input file [m] (scalar or table with 3 columns, h0x, h0y and h0)
%
% OUTPUT:
%    COAST
%      .h0     : Active profile height [m]
%      .h0_x   : Active profile x postition [m]
%      .h0_y   : Active profile y postition [m]
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2024 IHE Delft & Deltares
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

    if ischar(COAST.h0input)   % INPUT: from file
        % read table with active profile heights at point locations (with 3 columns, h0x, h0y and h0)
        % INPUT: EXAMPLE :   S.d='ActiveProfile.txt', 3 column text
        % file
        h0input=load(COAST.h0input);
        h0_x=COAST.h0input(:,1)';
        h0_y=COAST.h0input(:,2)';
        h0=COAST.h0input(:,3)';
        
        % find the right alongshore location for each of the active profile
        % interpolate h0         
        var1=struct;
        var2=struct;
        var2.h0_mc=h0;
        [~,var2i,~]=get_interpolation_on_grid('weighted_distance',COAST.x,COAST.y,h0_x,h0_y,var1,var2);        
        h0=var2i.h0_mc;  
        h0_x=COAST.x;
        h0_y=COAST.y;
        
    elseif isscalar(COAST.h0input)
        % Active profile at fixed height for the whole grid
        % INPUT: EXAMPLE :   COAST.PHIf0=10;
        h0_x=COAST.x;
        h0_y=COAST.y;
        h0=repmat(COAST.h0input,[1,COAST.n]);

    elseif isnumeric(COAST.h0input)  % array from input structure
        h0_x=COAST.h0input(:,1)';
        h0_y=COAST.h0input(:,2)';
        h0=COAST.h0input(:,3)';

        % find the right alongshore location for each of the active profile
        % interpolate h0         
        var1=struct;
        var2=struct;
        var2.h0_mc=h0;
        [~,var2i,~]=get_interpolation_on_grid('weighted_distance',COAST.x,COAST.y,h0_x,h0_y,var1,var2);        
        h0=var2i.h0_mc;  
        h0_x=COAST.x;
        h0_y=COAST.y;
        
    else     
        error('get_active_profile::Invalid active profile input');    
    end
    
    % OUTPUT: to structure
    COAST.h0=h0;
    COAST.h0_x=h0_x;
    COAST.h0_y=h0_y;
   
end
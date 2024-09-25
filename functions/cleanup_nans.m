function [COAST,i_diff]=cleanup_nans(COAST)
% function [COAST,i_diff]=cleanup_nans(COAST)
%
% cleans up the nans in the x,y coastline variables
%
% INPUT:
%    COAST
%      .x_mc        : x-coordinates for all coastal segments
%      .y_mc        : y-coordinates for all coastal segments
%
% OUTPUT:
%    COAST
%      .x_mc        : x-coordinates for all coastal segments (updated)
%      .y_mc        : y-coordinates for all coastal segments (updated)
%      .n_mc        : number of coastal segments (updated)
%    i_diff         : change in the number of points of x_mc/y_mc (-2 means two nan points were removed)
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

    i=1;
    eps=0.1;
    i_diff=0;
    while i<=length(COAST.x_mc)
        if i==1&& isnan(COAST.x_mc(i))    %[(Nan) 1 1 1 1] >> [1 1 1 1]
            COAST.x_mc=COAST.x_mc(2:end);
            COAST.y_mc=COAST.y_mc(2:end);
        elseif i==length(COAST.x_mc)&& isnan(COAST.x_mc(i))  %[1 1 1 (Nan)] >> [1 1 1 ]
            COAST.x_mc=COAST.x_mc(1:end-1);
            COAST.y_mc=COAST.y_mc(1:end-1);
        elseif i>1&&i<length(COAST.x_mc)&&isnan(COAST.x_mc(i-1))&&isnan(COAST.x_mc(i+1)) %[1 1 1 Nan (Nan) Nan 1 1] >> [1 1 1 1]
            COAST.x_mc=[COAST.x_mc(1:i-1),COAST.x_mc(i+2:end)];
            COAST.y_mc=[COAST.y_mc(1:i-1),COAST.y_mc(i+2:end)];
            if ~isnan(COAST.x_mc(i))
                i_diff=i_diff-1;  % administrate that 1 element is removed!
            end
        elseif i>1&&i<length(COAST.x_mc)-1&&isnan(COAST.x_mc(i-1))&&isnan(COAST.x_mc(i+2)) % [1 1 1 Nan (1) 1 Nan 1 1] >> [1 1 1 Nan 1 1]
            COAST.x_mc=[COAST.x_mc(1:i-1),COAST.x_mc(i+3:end)];
            COAST.y_mc=[COAST.y_mc(1:i-1),COAST.y_mc(i+3:end)];
            if ~isnan(COAST.x_mc(i)) || ~isnan(COAST.x_mc(i+1))
                i_diff=i_diff-1;                            % administrate that 1 element is removed!
            end
        elseif i==length(COAST.x_mc)-1&&isnan(COAST.x_mc(i))            % [1 1 (1)  Nan] >>  [1 1 1]
            COAST.x_mc=[COAST.x_mc(1:i-1)];
            COAST.y_mc=[COAST.y_mc(1:i-1)];
        elseif isnan(COAST.x_mc(i))&&isnan(COAST.x_mc(i+1))             % [1 1 1 (nan) nan 1 1 1] >>  %[1 1 1 (nan) 1 1 1]
            COAST.x_mc=[COAST.x_mc(1:i),COAST.x_mc(i+2:end)];
            COAST.y_mc=[COAST.y_mc(1:i),COAST.y_mc(i+2:end)];
        elseif i>1&&hypot(COAST.x_mc(i)-COAST.x_mc(i-1),COAST.y_mc(i)-COAST.y_mc(i-1))<eps  % remove repeated points    [.. xa xa ...] >> [.. xa ..] & [.. yb yb ..] >> [.. yb ..] 
            COAST.x_mc=[COAST.x_mc(1:i-1),COAST.x_mc(i+1:end)];
            COAST.y_mc=[COAST.y_mc(1:i-1),COAST.y_mc(i+1:end)];
            i=1;                                            % <- this means redo the whole loop to check if the remainder of the element is still valid & remove any remaining excessive nans -> possibly only an element with just one point is left, which needs to be removed by restarting the loop.           
        elseif i>2&&hypot(COAST.x_mc(i)-COAST.x_mc(i-2),COAST.y_mc(i)-COAST.y_mc(i-2))<eps ...  % if repeated points (same points with only one point in-between, possibly a triangular island, or local spike of coastline) [.. xa xb xa ...] >> [.. xa xb ..] & [.. yc yd yc ..] >> [.. yc yd ..]
                &&~isnan(COAST.x_mc(i-2))&&~isnan(COAST.x_mc(i-1))&&~isnan(COAST.x_mc(i))
            COAST.x_mc=[COAST.x_mc(1:i-1),COAST.x_mc(i+1:end)];
            COAST.y_mc=[COAST.y_mc(1:i-1),COAST.y_mc(i+1:end)];
            i=1;                                            % <- this means redo the whole loop to check if the remainder of the element is still valid & remove any remaining excessive nans -> possibly the 'repeated points' needs to be performed also in the restart of the loop to make it right again
        elseif i==2 && isnan(COAST.x_mc(i))                       % [1 Nan 1 1 1] >>  [1 1 1]  
            COAST.x_mc=[COAST.x_mc(i+1:end)];
            COAST.y_mc=[COAST.y_mc(i+1:end)];
            if ~isnan(COAST.x_mc(1))
                i_diff=i_diff-1;                            % administrate that 1 element is removed!
            end
        else
            i=i+1;
        end
    end
    
    COAST.n_mc=length(find(isnan(COAST.x_mc)))+1;
    
    %% make transport points xq_mc and yq_mc
    try
    [COAST]=get_transportpoints(COAST);
    end
end

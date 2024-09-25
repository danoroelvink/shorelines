function [props_mc,i1,i2]=insert_props(propsnew,props_mc,i_mc)
% function [props_mc,i1,i2]=insert_props(propsnew,props_mc,i_mc)
% 
% Insert updated properties props into multi-coastline properties
% Each row is a property, e.g. [x; y; s; Dfe; BW; Bf; Bm]
% props can be one matrix or a collection of row vectors, e.g. [x;y;s]
% [x_mc;y_mc;s_mc]=insert_props([xnew;ynew;snew],[x_mc;y_mc;s_mc],i_mc)
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

    nrows=size(propsnew,1);
    nans=zeros(nrows,1)*nan;
    %% Find indices i1 and i2, beginning and end of section in props_mc (only for first row)
    [ propsold,i1,i2,n_mc ] = get_one_polyvar( props_mc(1,:),i_mc );
    if i_mc<=n_mc
       %% Insert propsnew matrix (all rows) into props_mc
       props_mc=[props_mc(:,1:i1-1),propsnew,props_mc(:,i2+1:end)];
    else
       props_mc=[props_mc,nans,propsnew];
    end       
end

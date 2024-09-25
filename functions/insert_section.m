function [x_mc,y_mc,i1,i2]=insert_section(xnew,ynew,x_mc,y_mc,i_mc)
% function [x_mc,y_mc,i1,i2]=insert_section(xnew,ynew,x_mc,y_mc,i_mc)
%  OR
% function [x_mc,i1,i2]=insert_section(xnew,x_mc,i_mc)
% 
% This routine inserts sections of data (e.g. vectors of x-coordinates)
% into the overarching 'mc' variable, which contains data for all the 
% coastal elements (n_mc) at position of element i_mc.
% The existing data of element 'i_mc' is replaced by the new data. 
% It is possible to perform this action at once for two data variables
% with exactly the same size (e.g. x_mc and y_mc) or for just a single
% variable (e.g. QS_mc). 
%
% INPUT:
%     x_mc       : variable with a vector containing data of all the coastal segments separated by NaNs (e.g. x-coordinates)
%     y_mc       : (optional) additional variable with a vector containing other data of all the coastal segments separated by NaNs (e.g. y-coordinates) with exactly the same size as the x_mc
%     i_mc       : index of the element to be replaced in the x_mc / y_mc variables. 
%     xnew       : new data vector for element i_mc (e.g. y-coordinates for an element) with exactly the same size as the xnew
%     ynew       : (optional) additional new data vector for element i_mc (e.g. y-coordinates for an element) with exactly the same size as the xnew
% 
% OUTPUT: 
%     x_mc       : updated variable with a vector containing data of all the coastal segments separated by NaNs (e.g. x-coordinates)
%     y_mc       : (optional) additional updated variable with a vector containing other data of all the coastal segments separated by NaNs (e.g. y-coordinates) with exactly the same size as the x_mc
%     i_1        : index of the first point of the i_mc element in the original x_mc variable (up to this index-1 the original x_mc was preserved)
%     i_2        : index of the last point of the i_mc element in the original x_mc variable (from this index+1 the original x_mc was preserved)
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

    if nargin==3
        i_mc=x_mc;
        x_mc=ynew;
        y_mc=ynew;
        ynew=xnew;
    end

    %xnew(isnan(xnew))=-1e10;
    %ynew(isnan(ynew))=-1e10;
    [ xold,yold,n_mc,i1,i2 ] = get_one_polygon( x_mc,y_mc,i_mc );%% insert x,y back into COAST.x_mc,COAST.y_mc
    x_mc=[x_mc(1:i1-1),xnew,x_mc(i2+1:end)];
    y_mc=[y_mc(1:i1-1),ynew,y_mc(i2+1:end)];
    
    if nargout==3 && nargin==3
        y_mc=i1;
        i1=i2;
    end
    
end


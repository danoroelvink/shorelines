function [x_mc,y_mc,i1,i2]=insert_section(xnew,ynew,x_mc,y_mc,i_mc)
% function [x_mc,y_mc,i1,i2]=insert_section(xnew,ynew,x_mc,y_mc,i_mc)
% OR
% function [x_mc,i1,i2]=insert_section(xnew,x_mc,i_mc)
% 
% UNTITLED Summary of this function goes here
% Detailed explanation goes here
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


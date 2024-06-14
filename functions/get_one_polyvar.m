function [ v,i1,i2,n_mc] = get_one_polyvar( v_mc,i_mc )
% function [ v,i1,i2,n_mc] = get_one_polyvar( v_mc,i_mc )
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

    nans=find(isnan(v_mc));
    n_mc=length(nans)+1;
    if isempty(nans)
        i1=1;
        i2=length(v_mc);
    else
        if i_mc==1
            i1=1;
            i2=nans(i_mc)-1;
        elseif i_mc==n_mc;
            i1=nans(i_mc-1)+1;
            i2=length(v_mc);
        elseif i_mc>n_mc
            i1=0;
            i2=0;
            v=[];
            return;
        else
            i1=nans(i_mc-1)+1;
            i2=nans(i_mc)-1;
        end
    end
    v=v_mc(i1:i2);
    
end

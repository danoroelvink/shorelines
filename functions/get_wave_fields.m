function [ xg,yg,Hg,dirg,dirtab,Tptab ] = get_wave_fields_from_mat( fname )
% function [ xg,yg,Hg,dirg,dirtab,Tptab ] = get_wave_fields_from_mat( fname )
% 
% GET_WAVE_FIELDS
% routine to generate a series of wave fields for given wave
% height, wave directions and wave periods. The conditions are given  
% in jonswaptable.txt. Results computed by XBeach are stored in 
% xboutput.nc and read into the matrices xg,yg for the grid and
% Hg_all and phiwg_all; these represent a transformation matrix
% from which the wave field for an arbitrary wave condition can be
% interpolated.
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

    %% For now assume that the wave model has run and produced nHs by nphi 
    %% wave conditions
    load(fname)
    %xg=nc_varget(fname,'globalx');
    %yg=nc_varget(fname,'globaly');
    %Hg_all=nc_varget(fname,'H',[0 0 0],[Inf Inf Inf])*sqrt(2);
    %phiwg_all=nc_varget(fname,'thetamean',[0 0 0],[Inf Inf Inf]);

end

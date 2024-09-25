function [S]=initialize_randomgenerator(S)
% function [S]=initialize_randomgenerator(S)
%
% Set the seed of the random generator using S.randomseed
%
% INPUT:
%    S      : Settings of the ShorelineS model, with or without field '.randomseed'
%
% OUTPUT:
%    S      : Updated settings with field '.randomseedsettings'
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
    
    if (isoctave)
        if S.randomseed>=0
            % fixed seed
            S.randomseedsettings = rand('twister',get_twisterseed(S.randomseed));
        else
            % truely random
            S.randomseedsettings = rand('twister');
        end 
    else
        if S.randomseed>=0
            % fixed seed
            S.randomseedsettings = rng(S.randomseed);
        else
            % truely random
            S.randomseedsettings = rng('shuffle');
        end
    end
    
    %% write logfile
    struct2log(S,'S','w');
end
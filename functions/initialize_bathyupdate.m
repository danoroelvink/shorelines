function [BATHY]=initialize_bathyupdate(S)
% function [BATHY]=initialize_bathyupdate(S)
%
% INPUT:
%
% OUTPUT:
%
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

    fprintf('  Initialize bathy  \n');
    
    BATHY.bathy_update=S.bathy_update;
    BATHY.tide_interaction=S.tide_interaction;
    if ~isempty(S.bathy_update)  % For bathymetry update
        BATHY.update_time(:)=datenum(S.bathy_update(:,1),'yyyy-mm-dd');
        BATHY.tbu=1;
        BATHY.tupdate=BATHY.update_time(BATHY.tbu);
        BATHY.update_time(end+1)=0;
        %% write logfile
        % struct2log(BATHY,'BATHY','a');

    else
        BATHY.tupdate=[];
        BATHY.tbu=[];
        BATHY.update_time=[];
    end
end
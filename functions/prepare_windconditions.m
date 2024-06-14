function [WIND]=prepare_windconditions(S,TIME)
% function [WIND]=prepare_windconditions(S,TIME)
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

    WIND         = struct;
    WIND.dune    = S.dune;
    WIND.mud     = S.mud;
    WIND.dt      = S.dtdune;

    WND=[];
    if S.dune || S.mud
        WIND.rhoa    = S.rhoa;
        WIND.Cd      = S.Cd;
        WIND.uz      = S.uz;
        WIND.phiwnd  = S.phiwnd0;
        WIND.WNDfile = S.WNDfile;
        if isfield(S,'Windclimfile') && isempty(WIND.WNDfile)
            WIND.WNDfile = S.Windclimfile;
        end
        
        if ~isempty(WIND.WNDfile)
            if ischar(WIND.WNDfile)
                WIND.WNDfile={WIND.WNDfile};
            elseif size(WIND.WNDfile,2)>3 && (size(WIND.WNDfile,1)==1 || size(WIND.WNDfile,1)==3)
                WIND.WNDfile=WIND.WNDfile';
            end
            
            % read wind data files
            fprintf('  Read wind conditions for runup at dunes\n');
            [WND]=get_inputfiledata(WIND.WNDfile,TIME);
        end
        %% write logfile
        % struct2log(WIND,'WIND','a');
    end
    WIND.WND=WND;

end 

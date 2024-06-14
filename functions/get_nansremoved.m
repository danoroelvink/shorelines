function [x_mc,y_mc] = get_nansremoved(x_mc,y_mc)
% function [x_mc,y_mc] = get_nansremoved(x_mc,y_mc)
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

    nans=find(isnan(x_mc));
    x_mc(nans(diff(nans)==1))=-2e10;
    y_mc(nans(diff(nans)==1))=-2e10;
    x_mc=x_mc(x_mc~=-2e10);
    y_mc=y_mc(y_mc~=-2e10);
    if ~isempty(x_mc)
        if isnan(x_mc(1))
            x_mc=x_mc(2:end);
            y_mc=y_mc(2:end);
        end
    end
    if ~isempty(x_mc)
        if isnan(x_mc(end))
            x_mc=x_mc(1:end-1);
            y_mc=y_mc(1:end-1);
        end
    end

end


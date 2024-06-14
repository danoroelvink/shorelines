function [status]=runShorelineS(runfile);
% MODEL : ShorelineS
%
% This routine starts the ShorelineS model using an input file with all the settings.
% The model computes shoreline changes as a result of gradients in alongshore  
% sediment transport for arbitrary shaped coastlines. 
%
% INPUT:
%     runfile  filename of a textfile with datafields of ShorelineS input:
%              <fieldname1> = <property1>
%              <fieldname2> = <property2>
%              ...
%
% OUTPUT:
%     status   text output with the status of the ShorelineS run
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

    S=readkeys(runfile);
    [S,O] = ShorelineS(S);
    status=1;
end
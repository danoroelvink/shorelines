function [y]= get_randsample(n, k, wght)
% function [y]= get_randsample(n, k, wght)
%
% get_randsample Random wave condition using the probability of occurrence 
%
% INPUT:
%   n          number of values
%   k          should be 1
%   wght       probabilities per wave condition (in [%] or [days])
%
% OUTPUT:
%   y          randomly drawn sample using probability
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2022 IHE Delft & Deltares
%
%       Johan Reyns
%       j.reyns@un-ihe.org
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

    wght=wght(:)';
    if isempty(wght)
        y = randi(n,k,1);
    elseif length(wght)==n
        sumwght = sum(wght);
        p = wght / sumwght;
        pcum=min([0,cumsum(p)],1);  
        pcum(end) = 1;              % make sure the cumulative is exactly 1
        [~,y] = histc(rand(k,1),pcum);
    end
    y = y(:);

end

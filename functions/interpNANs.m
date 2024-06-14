function [data2]=interpNANs(data,useextrapolation)
% function [data2]=interpNANs(data)
%
% INPUT:
%     data                Array [Nx1] with data
%     useextrapolation    switch to use extrapolation (use 'nearest' or 'linear' for extrapolation)
%                         this does not affect the interpolation which will be 'linear' 
%
% data = [10, 9, nan, 8, nan, 7.5, nan, nan]
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

    if nargin<2
        useextrapolation='nearest';  %<- use extrapolation if it is set to 1
    end
    
    rotate=0;
    if max(size(data,1),2)<size(data,2)
        rotate=1;
        data=data';
    end
    
    if size(data,2)==1
        x = [1:length(data)]';
        y = data;
    else
        x = data(:,1);
        y = data(:,2);
    end
    
    % make sure that x is increasing
    [x,ids]  = sort(x);
    y        = y(ids);
    eps      = 1e-6;
    for tt=1:length(x)-1
        if x(tt)==x(tt+1)
            x(tt+1)=x(tt+1)+eps;
        end
    end
    
    %% remove nans
    id = find(~isnan(y));
    x2 = x(id);
    y2 = y(id);
    
    if length(id)==length(x)
         ynew = y2;
    else
        if length(x2)>1
            if strcmpi(useextrapolation,'linear')
                ynew = interp1(x2,y2,x,'linear','extrap');
            else
                y3 = interp1(x2,y2,x,'linear');
                id = find(~isnan(y3));
                x2 = x(id);
                y2 = y3(id);
                ynew = interp1(x2,y2,x,'nearest','extrap');
            end
        elseif length(x2)==1
            ynew=repmat(y2,size(x));
        else
            %fprintf('Need more data points than 1\n')
            ynew=nan(size(x));
        end
    end
    
    if size(data,2)==2
        data2 = [x,ynew];
    else
        data2 = [ynew];
    end
    
    if rotate==1
        data2 = data2';
    end
end

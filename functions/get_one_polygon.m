function [ x,y,n_mc,i1,i2 ] = get_one_polygon( x_mc,y_mc,i_mc )
% function [ x,n_mc,i1,i2 ] = get_one_polygon( x_mc,i_mc )
% function [ x,y,n_mc,i1,i2 ] = get_one_polygon( x_mc,y_mc,i_mc )
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

    if nargin<3    
        i_mc=y_mc;
        y_mc=x_mc;
    end
    
    nans=find(isnan(x_mc));
    nansremove=nans([false,diff(nans)==1]);
    %nansremove=nans(diff(nans)==1); 
    
    if ~isempty(nansremove)
        idkeep=setdiff([1:length(x_mc)],nansremove);
        x_mc=x_mc(idkeep);
        y_mc=y_mc(idkeep);
        nans=find(isnan(x_mc));
    end   
    
    n_mc=length(nans)+1;
    %i_mc=min(i_mc,n_mc);
    if isempty(nans)
        i1=1;
        i2=length(x_mc);
    else
        if i_mc==1
            i1=1;
            i2=nans(i_mc)-1;
        elseif i_mc==n_mc
            i1=nans(i_mc-1)+1;
            i2=length(x_mc);
        elseif i_mc>n_mc
            %disp('i_mc>n_mc')
            i1=nans(n_mc-1)+1;
            i2=length(x_mc);
        else
            i1=nans(i_mc-1)+1;
            i2=nans(i_mc)-1;
        end
    end
    x=x_mc(i1:i2);
    y=y_mc(i1:i2);
        
    if isempty(x_mc(~isnan(x_mc)))
        n_mc=0;
    end
    
    if nargin<3 && nargout==4   
        y=n_mc;
        n_mc=i1;
        i1=i2;
    end

end

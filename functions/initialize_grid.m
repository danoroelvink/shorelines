function [x_mc,y_mc]=initialize_grid(x_mc,y_mc,ds0)
% function [x_mc,y_mc]=initialize_grid(x_mc,y_mc,ds0)
% 
% find number of sections
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
    n_mc=length(nans)+1;
    yesplot=false;
    if yesplot
       figure(11)
    end
    for i_mc=1:n_mc
        [ x,y,n_mc,i1,i2 ] = get_one_polygon( x_mc,y_mc,i_mc );   
        if yesplot
           plot(x,y); hold on
        end
        % remove double points
        IDunique = ~[0,(diff(x)==0&diff(y)==0)];
        x0=x(IDunique);
        y0=y(IDunique);
        % compute distance along line
        s0=zeros(size(x0));
        for i=2:length(x0)
            s0(i)=s0(i-1)+hypot(x0(i)-x0(i-1),y0(i)-y0(i-1));
        end
        ns=ceil(s0(end)/ds0);
        ds=s0(end)/ns;
        s=[0:ds:s0(end)];
        % interpolate x-values
        if length(s)>1
            x=interp1(s0,x0,s);
            y=interp1(s0,y0,s);
            %% insert x,y back into x_mc,y_mc
            x_mc=[x_mc(1:i1-1),x,x_mc(i2+1:end)];
            y_mc=[y_mc(1:i1-1),y,y_mc(i2+1:end)];
        else
            n_mc=1;
        end
    end
    if yesplot
       plot(x_mc,y_mc,'o')
    end
end

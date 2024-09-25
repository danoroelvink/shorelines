function [xi,yi]=select_multi_polygon(col)
% function [xi,yi]=select_multi_polygon(col)
% 
% Select polygon to include in bathy
%
% INPUT:
%     col         : colour used for the plot (e.g. 'r')
%
% OUTPUT:
%     xi          : x-coordinates of selected polygon [m]
%     yi          : y-coordinates of selected polygon [m]
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

    hold on;
    xi = [];yi=[];
    n = 0;
    % Loop, picking up the points.
    finished=0;
    while ~finished
        nold=n;
        but = 1;
        while but == 1
            [xs,ys,but] = ginput(1);
            if but==1
                n = n+1;
                xi(n)=xs;
                yi(n)=ys;
                plot(xi,yi,[col '-o']);
            end
        end
        if nold==n || but==113 || but==27                                      % respectively 'q' or 'escape'
            finished=1;
        else
            n=n+1;
            xi(n)=nan;
            yi(n)=nan;
        end
    end
    
    % remove single points!
    idnan = find(isnan(xi));
    xi2=[];
    yi2=[];
    if length(idnan)>=1
        if idnan(1)>1
            idnan=[0;idnan(:)]';
        end
        if ~isnan(xi(n))
            xi=[xi,nan];
            yi=[yi,nan];
            n=n+1;
            idnan=[idnan,n];
        end
        id = diff(idnan)-1;
        for mm=1:length(id)
            if id(mm)>1
                xi2 = [xi2,xi(idnan(mm)+1:idnan(mm+1))];
                yi2 = [yi2,yi(idnan(mm)+1:idnan(mm+1))];
            end
        end
        xi=xi2;
        yi=yi2;
        n=length(xi);
    end
    
    % remove trailing nan
    if n>0&&isnan(xi(n))
        xi=xi(1:n-1);
        yi=yi(1:n-1);
    end
    
end

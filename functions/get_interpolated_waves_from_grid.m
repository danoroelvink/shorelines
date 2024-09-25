function [H,phiw_cd,x_cd,y_cd]=get_interpolated_waves_from_grid(x,y,n,nq,surfwidth,xg,yg,Hg,phiwg)
% function [H,phiw_cd,x_cd,y_cd]=get_interpolated_waves_from_grid(x,y,n,nq,surfwidth,xg,yg,Hg,phiwg)
% 
% This function interpolates waves from a grid.
% 
% INPUT: 
%    x              : x-coordinates of the coastline [m]
%    y              : y-coordinates of the coastline [m]
%    n              : number of grid cells of the coastline
%    nq             : number of grid cells of the qs-points
%    surfwidth      : surfzone width [m]
%    xg             : x-coordinates of the wave grid [m]
%    yg             : y-coordinates of the wave grid [m]
%    Hg             : wave height field of the wave grid [m]
%    phiwg          : wave direction fields of the wave grid [°N]
%
% OUTPUT: 
%    H              : wave height [m]
%    phiw_cd        : wave direction [°N]
%    x_cd           : x-coordinates of wave stations [m]
%    y_cd           : y-coordinates of wave stations [m]
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

    H=[];
    phiw_cd=[];
    %% create line at depth of closure x_cd, y_cd
    %COAST.n=length(COAST.x)-1;
    for i=1:nq
        i1=mod(i-1,n);
        i2=mod(i,n);
        
        dX=x(i2)-x(i1);
        dY=y(i2)-y(i1);
        Hyp=hypot(dX,dY);
        dx=surfwidth*sind(phiwg(1,1));
        dy=surfwidth*cosd(phiwg(1,1));
        x_cd(i)=0.5*(x(i1)+x(i2))+dx;
        y_cd(i)=0.5*(y(i1)+y(i2))+dy;
    end
    for i=1:nq
        if ~isnan(x)
            [row,col]=find_in_grid(xg,yg,x_cd(i),y_cd(i));
            if (size(row,1)>1)
                break %continue%/break to check
            end
            H(i)=Hg(row,col);
            phiw_cd(i)=phiwg(row,col);  
        end
    end
    phiw_cd(isnan(phiw_cd))=phiwg(1,1);
end
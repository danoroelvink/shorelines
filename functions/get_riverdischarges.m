function [MUD]=get_riverdischarges(TIME,COAST,TRANSP,MUD)
% function [MUD]=get_riverdischarges(TIME,COAST,TRANSP,MUD)
%
% INPUT: 
%    TIME
%       .tnow          : current date in the model (in timenum format)
%    COAST
%       .x             : x-coordinates of the coastline [m]
%       .y             : y-coordinates of the coastline [m]
%       .n             : Number of coastline points (length of x and y)
%    MUD
%       .tstart        : Start of riverdischarges
%       .tend          : End of riverdischarge 
%       .x             : x-coordinates of riverdischarges
%       .y             : y-coordinates of riverdischarges
%       .n             : Number of riverdischarges
%       .rate          : riverdischarge rate for each individual measure in [m3/yr]
%
% OUTPUT:
%     MUD              :     
%       .rate_density  : actual nourishment rate in [m3/m/year]
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

    MUD.rate_density=zeros(1,COAST.n);
    if MUD.used
        for i_n=1:MUD.nriv
            if TIME.tnow>=MUD.tstart(i_n) && TIME.tnow<=MUD.tend(i_n)
                % find indices of begin and end point of each river mouth
                % compute number of grid cells in-between begin and end point
                dxmud=[];
                for i_mc=1:COAST.n_mc
                    [ xc,yc,~,~,~ ] = get_one_polygon( COAST.x_mc,COAST.y_mc,i_mc );
                    dist1=hypot(xc-MUD.xriv(i_n,1),yc-MUD.yriv(i_n,1));
                    dist2=hypot(xc-MUD.xriv(i_n,2),yc-MUD.yriv(i_n,2));
                    idx1=find(dist1==min(dist1),1);
                    idx2=find(dist2==min(dist2),1);
                    idx{i_mc}=[min(idx1,idx2):max(idx1,idx2)];
                    % Check if coastal section is cyclic
                    cyclic=get_cyclic(xc,yc,COAST.ds0);
                    % Compute grid cell size
                    s=[0,cumsum((diff(xc).^2+diff(yc).^2).^0.5)];
                    ds0=diff(s);
                    ds=[ds0(1),(ds0(1:end-1)+ds0(2:end))/2,ds0(end)]; 
                    if COAST.cyclic(i_mc)
                        ds(1)=(ds0(1)+ds0(end))/2;
                        ds(end)=(ds0(1)+ds0(end))/2; 
                    end
                    dxmin=1e-5;
                    dxmud(i_mc)=max(sum(ds(idx{i_mc})),dxmin); 
                    if isnan(dxmud(i_mc))
                        dxmud(i_mc)=0;
                    end
                end 
                
                % add river discharge
                % mud [m3/m/year]
                % make sure to use only the fraction that is nourished at the considered coastal segment i_mc
                fraction = dxmud./sum(dxmud);
                if ~isempty(idx{COAST.i_mc}) 
                    MUD.rate_density(idx{COAST.i_mc})=MUD.rate_density(idx{COAST.i_mc})+MUD.rate(i_n)./dxmud(COAST.i_mc).*fraction(COAST.i_mc);
                end
            end
        end
    end
end  

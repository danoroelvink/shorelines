function [fraction,i_mc_rel,idx,dxnour] = get_fnourishment_fraction(xn,yn,x_mc,y_mc,n_mc,ds_mc,shadow_mc)
% function [fraction,i_mc_rel,idx,dxnour] = get_fnourishment_fraction(xn,yn,x_mc,y_mc,n_mc,ds_mc,shadow_mc)
%
% Computes the fraction of the total nourishment volume per coastal segment
%
% INPUT: 
%   xn
%   yn
%   x_mc
%   y_mc
%   n_mc
%   ds_mc
%   shadow_mc
%
% OUTPUT: 
%   fraction : fraction of total nourishment volume per coastal segment [-]
%   i_mc_rel : relevant segments for nourishment
%   idx      : indices of coastline affected by nourishment
%   dxnour   : nourishment length per segment [m]
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

    dxnour   = [];
    L1       = [];
    L2       = [];
    fraction = [];
                
    for i_mc = 1 : n_mc % loop over nourishments in main function 

        xc             = get_one_polygon(x_mc,i_mc);
        yc             = get_one_polygon(y_mc,i_mc);
        dsc            = get_one_polygon(ds_mc,i_mc);
        dist1          = ((xc-xn(1)).^2+(yc-yn(1)).^2).^0.5;
        dist2          = ((xc-xn(2)).^2+(yc-yn(2)).^2).^0.5;
        idx1           = find(dist1==min(dist1),1);
        idx2           = find(dist2==min(dist2),1);
        idx{i_mc}      = [min(idx1,idx2):max(idx1,idx2)];
                    
        % Check if coastal section is cyclic
        dsc(isnan(dsc))=0;
        ds0=median(dsc);
        cyclic = get_cyclic(xc,yc,ds0);
        if cyclic 
            if idx1 > idx2
                idx{i_mc} = [ 1 : idx2  idx1 : length(dist1) ] ;
            end 
        end 

        % Check for shadow zones
        [ shadowS ]    = get_one_polygon(shadow_mc,i_mc);
        id_shadow      = find( shadowS == 1 ); 
        Lia            = ismember( idx{i_mc} , id_shadow );
        idx{i_mc}(Lia) = []; 
        
        L1(i_mc)  = min(dist1);
        L2(i_mc)  = min(dist2);

        dxnour(i_mc) = sum(dsc(idx{i_mc})); 
        if isnan(dxnour(i_mc))
           dxnour(i_mc) = 0;
        end
    
    end 
    
    % Find relevant segments
    [~,id1] = min(L1);
    [~,id2] = min(L2);
    i_mc_rel  = [id1 id2];
    if length(unique(i_mc_rel)) == 1 % fix for one segment 
        i_mc_rel = i_mc_rel(1); 
    end 
    i_mc_del  = 1 : n_mc; 
    [val,ia,ib] = intersect(i_mc_rel,i_mc_del);
    i_mc_del(ib) = [];
    idx(i_mc_del) = {[]}; 
    dxnour(i_mc_del) = 0; 
    
    % Compute fraction 
    fraction = dxnour./sum(dxnour);
    
end 

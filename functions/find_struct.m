function [STRUC,TRANSP]=find_struct(COAST,STRUC,TRANSP)
% function [STRUC,TRANSP]=find_struct(COAST,STRUC,TRANSP);
%
% Finds the indices of structures on the coastline (for xy-points)
%
% INPUT:
%   COAST
%      .x       : x-coordinates of coastal section
%      .y       : y-coordinates of coastal section
%   STRUC
%     .x_hard   : x-points of hard structures (multiple structures are seperated with NAN's)
%     .y_hard   : y-points of hard structures (multiple structures are seperated with NAN's)
%
% OUTPUT:
%   STRUC.structS     : Indices of structures on the coastline (for xy-points)
%   TRANSP.QS         : Transport set to 0 at 'strucS' indices
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


    x=COAST.x;
    y=COAST.y;    
    structS=zeros(1,length(x)-1);
    structS=logical(structS);
    
    if ~isempty(STRUC.x_hard)
        eps=1e-6;
        [xx,yy,indc,inds]=get_intersections(x,y,STRUC.x_hard,STRUC.y_hard);
        structS(floor(indc+eps))=true;
        if length(indc)>2
            %figure;plot(x,y,'k.-');hold on;plot(STRUC.x_hard,STRUC.y_hard,'r.-')
        end
        if isfield(COAST,'x')
            structS(1)=0;
            structS(end)=0; % avoid problems with real groynes
            TRANSP.QS(structS)=0;
        end
    end    
    
    STRUC.structS=structS;
end

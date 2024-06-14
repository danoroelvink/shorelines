function [DELTA]=prepare_delta(S)
% [DELTA]=prepare_delta(S)
%
% INPUT:
%    S
%         .flood_delta    :  
%         .x_flood_pol    :  
%         .y_flood_pol    :  
%         .x_spit_pol     :  
%         .y_spit_pol     :  
%         .dxf            :  
%         .overdepth      :  
%         
% OUTPUT:
%    DELTA
%         .x_flood
%         .y_flood
%         .flood_deficit
%         .fcell_area
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

    fprintf('  Prepare delta \n');
    
    %% Flood delta
    DELTA=struct;
    DELTA.x_flood=[];
    DELTA.y_flood=[];
    DELTA.flood_deficit=[];
    DELTA.fcell_area=[];
    
    if S.flood_delta & S.channel
        if ~isempty(S.x_flood_pol)
            x_flood_pol=S.x_flood_pol;
            y_flood_pol=S.y_flood_pol;
            x_spit_pol=S.x_spit_pol;
            y_spit_pol=S.y_spit_pol;
        elseif ~isempty(S.LDBflood)
            xy_flood=load(S.LDBflood);
            x_flood_pol=xy_flood(:,1)'-S.XYoffset(1);
            y_flood_pol=xy_flood(:,2)'-S.XYoffset(2);
            xy_spit=load(S.LDBspit);
            x_spit_pol=xy_spit(:,1)'-S.XYoffset(1);
            y_spit_pol=xy_spit(:,2)'-S.XYoffset(2);
        else
            figure;plot(x_mc,y_mc,xr_mc,yr_mc,'--');axis equal;hold on;
            disp('select flood delta outline');
            [x_flood_pol,y_flood_pol]=select_multi_polygon('k');
            disp('select spit outline');
            [x_spit_pol,y_spit_pol]=select_multi_polygon('k');
            save('flood_delta.mat','x_flood_pol','y_flood_pol','x_spit_pol','y_spit_pol');
        end
        [DELTA.x_flood,DELTA.y_flood,DELTA.flood_deficit,DELTA.fcell_area] = prepare_flood_delta(x_flood_pol,y_flood_pol,S.dxf,S.dxf,S.overdepth);
        %% write logfile
        % struct2log(DELTA,'DELTA','a');

    end
end

function [DELTA]=prepare_delta(S)
% function [DELTA]=prepare_delta(S)
%
% Prepares the DELTA data-structure.
%
% INPUT: 
%    S
%         .flooddelta    : switch (0/1) for flood delta 
%         .ldbflood      : file with x and y coordinates of the flood delta polygon [Nx2] or directly provide a matrix
%         .xfloodpol     : x-coordinates of flood delta polygon
%         .yfloodpol     : y-coordinates of flood delta polygon
%         .ldbspit       : file with x and y coordinates of the spit polygon [Nx2] or directly provide a matrix
%         .xspitpol      : x-coordinates of the spit polygon
%         .yspitpol      : y-coordinates of the spit polygon
%         .dxf           : resolution of flood delta area [m]
%         .overdepth     : initial overdepth flood delta [m]
%         
% OUTPUT:
%    DELTA
%         .xflood        : x-coordinates of flood delta polygon
%         .yflood        : y-coordinates of flood delta polygon
%         .flooddeficit  : flood delta volume deficit
%         .fcellarea     : flood delta area
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

    fprintf('  Prepare delta \n');
    
    %% Flood delta
    DELTA=struct;
    DELTA.xflood=[];
    DELTA.yflood=[];
    DELTA.flooddeficit=[];
    DELTA.fcellarea=[];
    
    if S.flooddelta & S.channel
        if ~isempty(S.xfloodpol)
            xfloodpol=S.xfloodpol;
            yfloodpol=S.yfloodpol;
            xspitpol=S.xspitpol;
            yspitpol=S.yspitpol;
        elseif ~isempty(S.ldbflood)
            xyflood=load(S.ldbflood);
            xfloodpol=xyflood(:,1)'-S.xyoffset(1);
            yfloodpol=xyflood(:,2)'-S.xyoffset(2);
            xyspit=load(S.ldbspit);
            xspitpol=xyspit(:,1)'-S.xyoffset(1);
            yspitpol=xyspit(:,2)'-S.xyoffset(2);
        else
            figure;plot(x_mc,y_mc,xrmc,yrmc,'--');axis equal;hold on;
            disp('select flood delta outline');
            [xfloodpol,yfloodpol]=select_multi_polygon('k');
            disp('select spit outline');
            [xspitpol,yspitpol]=select_multi_polygon('k');
            save('flooddelta.mat','xfloodpol','yfloodpol','xspitpol','yspitpol');
        end
        [DELTA.xflood,DELTA.yflood,DELTA.flooddeficit,DELTA.fcellarea] = prepare_flooddelta(xfloodpol,yfloodpol,S.dxf,S.dxf,S.overdepth);
        
        %% write logfile
        % struct2log(DELTA,'DELTA','a');
    end
end

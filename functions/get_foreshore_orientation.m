function [COAST]=get_foreshore_orientation(COAST)
% function [COAST]=get_foreshore_orientation(COAST)
% 
% INPUT: 
%   COAST        : Structure with data of the ShorelineS model, of which is used
%     .phif0     : user input for the shoreline orientation [°]
%     .x         : x-coordinate of coastal segment [m]
%     .y         : y-coordinate of coastal segment [m]
%     .n         : number of grid cells of coastal segment
%     .PHIc      : Coastline orientation at each grid cell [°]
%
% OUTPUT:
%   COAST
%     .PHIf      : Lower shoreface orientation at each grid cell for the current coastal segment [°]
%     .PHIf_x    : x-coordinate at each grid cell for all of the coastal segments [m]
%     .PHIf_y    : y-coordinate at each grid cell for all of the coastal segments [m]
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2021 IHE Delft & Deltares
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
            
    % shoreface orientation (PHIf)    
    if ischar(COAST.PHIf0)   % INPUT: from file
        % read table with predefined shoreface orientation (with 3 columns, PHIx, PHIy and PHIf orientation)
        % INPUT: EXAMPLE :   S.phi='foreshoreorientation.txt', 3 column text
        % file
        COAST.PHIf0=load(COAST.PHIf0);
    end
    
    if isempty(COAST.PHIf0) || size(COAST.PHIf0,2)==4 % INPUT: at t=0
        % lower shoreface at original orientation of coastline at t=0
        % NO INPUT      : use coastline orientation (e.g. in case of CERC)
        % this is only done at t=0. Afterwards the assessed PHIf at t=0 is reinterpolated on the grid every time step (i.e. when field 'PHIf' is available). 
        
        % do a little bit of smoothing with just the next cell
        nl=length(COAST.PHIc);
        xq=COAST.xq(:);
        yq=COAST.yq(:);
        if COAST.cyclic
            %PHIfsmooth=COAST.PHIc(2:end-1);
            PHIf=get_smoothdata(repmat(COAST.PHIc(2:end-1),[1 3]),'angle',1);
            PHIfsmooth=mod(PHIf(nl-1:2*nl-4),360);   
            idx=[2:nl-1];
            PHIfsmooth=mod(PHIf(nl-1:2*nl-3),360);   
            idx=[2:nl];
        else
            PHIf=get_smoothdata(COAST.PHIc,'angle',1);           
            PHIfsmooth=mod(PHIf,360);
            idx=[1:nl];
        end
        
        % add foreshore orientations to variable until foreshore orientation is determined for all coast sections
        if COAST.i_mc==1
            COAST.PHIf0=[xq(idx),yq(idx),PHIfsmooth(:),PHIfsmooth(:)];
        else
            COAST.PHIf0=[COAST.PHIf0;[xq(idx),yq(idx),PHIfsmooth(:),PHIfsmooth(:)]];
        end
        if COAST.i_mc==COAST.n_mc
            COAST.PHIf0=COAST.PHIf0(:,1:3);
        end
    end
    
    if iscell(COAST.PHIf0) && length(COAST.PHIf0)>=2 % INPUT: at t=0
        % lower shoreface angles based on smoothed coastline (only present section in COAST.PHIc)
        % INPUT: EXAMPLE :   COAST.PHIf0={'gaussian',7};
        smoothmethod=COAST.PHIf0{1};
        smoothrange=COAST.PHIf0{2};
        if COAST.cyclic
            nl=length(COAST.PHIc);
            cPHIf=smoothdata(cosd(repmat(COAST.PHIc,[1 3])),smoothmethod,smoothrange);
            sPHIf=smoothdata(sind(repmat(COAST.PHIc,[1 3])),smoothmethod,smoothrange);
            PHIf=atan2d(sPHIf(nl+1:2*nl),cPHIf(nl+1:2*nl));
        else
            cPHIf=smoothdata(cosd(COAST.PHIc),smoothmethod,smoothrange);     
            sPHIf=smoothdata(sind(COAST.PHIc),smoothmethod,smoothrange);
            PHIf=atan2d(sPHIf,cPHIf);
        end
        PHIfsmooth=mod(PHIf,360);
        COAST.PHIf=PHIfsmooth;
        
        if COAST.i_mc==1
            COAST.PHIf0={COAST.PHIf0{1},COAST.PHIf0{2},COAST.xq(:),COAST.yq(:),PHIfsmooth(:)};
        else
            COAST.PHIf0{3}=[COAST.PHIf0{3};COAST.xq(:)];
            COAST.PHIf0{4}=[COAST.PHIf0{4};COAST.yq(:)];
            COAST.PHIf0{5}=[COAST.PHIf0{5};PHIfsmooth(:)];
        end
        if COAST.i_mc==COAST.n_mc
            COAST.PHIf0=[COAST.PHIf0{3},COAST.PHIf0{4},COAST.PHIf0{5}];
        end
    end
    
    if isscalar(COAST.PHIf0)
        % lower shoreface at fixed orientation for the whole grid
        % INPUT: EXAMPLE :   COAST.PHIf0=312;
        COAST.PHIf=repmat(COAST.PHIf0,[1,length(COAST.xq)]);
    end
    
    % lower shoreface at predefined orientation in table
    % INPUT: EXAMPLE :   PHIf0=[x1,y1,phif1; ... ; xn,yn,phifn]       
    if ~isscalar(COAST.PHIf0) && ~iscell(COAST.PHIf0)
        PHIf_x=COAST.PHIf0(:,1)';
        PHIf_y=COAST.PHIf0(:,2)';
        PHIf0=COAST.PHIf0(:,3)';
        
        % find the right alongshore location for each of the wave climates
        % re-interpolate PHIf          
        if COAST.cyclic
            [~,PHIf1,~]=get_interpolation_on_grid('weighted_distance',COAST.xq,COAST.yq,PHIf_x,PHIf_y,[],PHIf0);
        else
            [~,PHIf1,~]=get_interpolation_on_grid('alongshore_mapping',COAST.xq,COAST.yq,PHIf_x,PHIf_y,[],PHIf0);
        end
        COAST.PHIf=PHIf1;
    end
    
end

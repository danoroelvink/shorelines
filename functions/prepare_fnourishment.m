function [FNOUR] = prepare_fnourishment(S);
% function [FNOUR] = prepare_fnourishment(S)
%
% Prepares all required input for the computation of the shoreface nourishment source term.
%
% INPUT:  
%   S
%       .fnourish : switch for using shoreface nourishments (0/1)
%       .mb       : coefficient 1 for shoreface nourishments supply over time
%       .labda0   : coefficient 2 for shoreface nourishments supply over time
%       .rhos     : density of the sediment [kg/m3]
%       .sal      : salinity of the water [ppm salt]
%       .temp     : temperature of the water [degrees celsius]
%       .fnorfile : shoreface nourishment input file 
% 
% OUTPUT:
%   FNOUR
%       .fnourish : switch for using shoreface nourishments (0/1)
%       .x        : x-coordinates of the start-end points of the shoreface nourishments [m]
%       .y        : y-coordinates of the start-end points of the shoreface nourishments [m]
%       .L        : length of the shoreface nourishments [m]
%       .d        : depth of the shoreface nourishments [m]
%       .d50      : median grain size of the shoreface nourishments [m]
%       .tstart   : moment of implementation of the shoreface nourishments [days in datenum format]
%       .V        : volume of the shoreface nourishments
%       .n        : number of shoreface nourishments
%       .mb       : coefficient 1 for shoreface nourishments supply over time
%       .labda0	  : coefficient 2 for shoreface nourishments supply over time
%       .rhos     : density of the sediment [kg/m3]
%       .sal      : salinity of the water [ppm salt]
%       .temp     : temperature of the water [degrees celsius]
%       .K        : diffusion coefficients for the shoreface nourishments, on the basis of alongshore volumetric density of the nourishments
%       .w        : fall velocities for each of the shoreface nourishments [m/s]
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

fprintf('  Prepare shoreface nourishments \n');

% Initialize struct 
FNOUR           = struct;
FNOUR.x         = [];
FNOUR.y         = [];
FNOUR.L         = [];
FNOUR.d         = [];
FNOUR.d50       = [];
FNOUR.tstart    = [];
FNOUR.V         = [];
FNOUR.n         = [];
FNOUR.K         = []; 
FNOUR.V_t       = []; 
FNOUR.fnourish  = S.fnourish; % check if shoreface nourishments are activated

% Get info from default values 
FNOUR.mb     = S.mb; 
FNOUR.labda0 = S.labda0;
FNOUR.rhos   = S.rhos;
FNOUR.sal    = S.sal; 
FNOUR.temp   = S.temp; 

if ~isempty(FNOUR.fnourish) && FNOUR.fnourish ~= 0 % if shoreface nourishments are activated

    if ~isempty(findstr(lower(S.fnorfile),'.fnor')) % check if fnorfile is present 
        
        % fnorfile contains [xstart, xend, ystart,  yend, L, d, d50, tstart(yyyymmdd), totalvolume]
        fnor            = load(S.fnorfile);
        FNOUR.x         = [fnor(:,1) fnor(:,2)];
        FNOUR.y         = [fnor(:,3) fnor(:,4)];
        FNOUR.L         = fnor(:,5);
        FNOUR.d         = fnor(:,6); 
        FNOUR.d50       = fnor(:,7); 
        FNOUR.tstart    = datenum(num2str(fnor(:,8)),'yyyymmdd');
        FNOUR.V         = fnor(:,9); % [m3] 
        FNOUR.n         = size(fnor,1);
 
        for i = 1:FNOUR.n

            % Compute K 
            if isempty(S.k)
                FNOUR.K(i) = get_fnourishment_diffusion( FNOUR.V(i) / FNOUR.L(i) ); % K-formulation derived for density  
            else 
                FNOUR.K(i) = S.k(i);
            end 
            
            % Compute fall velocity 
            FNOUR.w(i) = get_fallvelocity(FNOUR.d50(i), FNOUR.rhos, FNOUR.sal, FNOUR.temp);
        end 
       
    else
        fprintf('  No shoreface nourishment info file present \n');
    end
    %% write logfile
    % struct2log(FNOUR,'FNOUR','a');

end
 












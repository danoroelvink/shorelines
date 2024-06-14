function [Vend,Qend,Fac] = get_fnourishment_rate(dt,dtsteps, V0, labda0, w, mb, K, WAVE, idx_mc, x_fnour, y_fnour)
% function [Vend,Qend,Fac] = get_fnourishment_rate(dt,dtsteps, V0, labda0, w, mb, K, WAVE, idx_mc, x_fnour, y_fnour)
%
% Get the relative volume of a shoreface nourishment still present at the
% starting time of the simulation if starting time nourishment < starting
% time simulation
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

    %% Get average wave conditions over the considered 1) time period and 2) spatial impact area
    if dtsteps==1 
        % default method used for evaluations per timestep after model start time (t0)
        H = sqrt( mean( WAVE.HSo(idx_mc).^2 ) );
        T = sqrt( mean( WAVE.TP(idx_mc).^2 ) ); 

    else
        % for evaluation of the conditions before t0 (i.e. when initializing V0)
        if isfield(WAVE.WVC,'Hs') && size(WAVE.WVC,2) == 1 && dtsteps>1 
            % timeseries, alongshore uniform 
            H = sqrt( mean( WAVE.WVC.Hs.^2 ) ); 
            T = sqrt( mean( WAVE.WVC.Tp.^2 ) ); 

        elseif isfield(WAVE.WVC,'Hs') && size(WAVE.WVC,2) > 1 && dtsteps>1 
            % timeseries, alongshore varying 
            x_wvc = [WAVE.WVC(:).x]; 
            y_wvc = [WAVE.WVC(:).y];
            x_fnour = mean( x_fnour ); % middle of nourishment 
            y_fnour = mean( y_fnour );             
            dist = sqrt( ( x_wvc - x_fnour ).^2 + ( y_wvc - y_fnour ).^2 ); % Find closest point to nourishment 
            [~,id] = min(dist);
            H = sqrt( mean( WAVE.WVC(id).Hs.^2 ) ); 
            T = sqrt( mean( WAVE.WVC(id).Tp.^2 ) ); 

        elseif isempty(WAVE.WVC)
            % uniform conditions in time and space
            H = sqrt( mean( WAVE.HSo(idx_mc).^2 ) );
            T = sqrt( mean( WAVE.TP(idx_mc).^2 ) ); 

        end
    end 

    % Compute volume change of the shoreface nourishment over time
    V(1) = V0; % volume density of the foreshore nourishment [m3/m]
    q(1) = 0;
    for t = 1:max(dtsteps,1)
        t2=t+1;
        labda    = labda0 * ( H / (w*T ) ) ^ mb;            % response coefft2cient [-]
        q(t2)    = ( V(t) * (1-exp(-1*K*labda*dt)) ) / dt;  % rate of supply of the nourishment to the beach [m3/m/s]
        dVdt(t2) = -1*q(t2);                                % rate of erosion of the nourishment [m3/m/s]
        V(t2)    = V(t) + dVdt(t2)*dt;                      % new volume density of the foreshore nourishment [m3/m]
    end 

    % Get factor 
    Qend = q(end);        % [m3/m/s]
    Vend = V(dtsteps+1);  % [m3/m]
    Fac  = Vend / V(1);   % [-]

end



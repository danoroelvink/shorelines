function [TIME]=get_timestep(TIME,COAST,TRANSP)
% function [adt]=get_timestep(adt,dt,automatic,s,d,QS,QSmax)
% 
% INPUT:
%   TIME
%      .adt
%      .dt
%   COAST
%      .s
%      .d
%      .h0
%   TRANSP
%      .QS
%      .QSmax 
%
% OUTPUT:
%   adt
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
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses
%   --------------------------------------------------------------------

    if ~TIME.tc
        adt=TIME.dt;
    else
        dsmin=min(diff(COAST.s));            
        h0min=min(COAST.h0);
        adt=dsmin^2*h0min / (4*max(abs(TRANSP.QSmax)));
        if max(abs(TRANSP.QS))*4<max(abs(TRANSP.QSmax))
            adt=0.25*dsmin^2*h0min / (4*max(abs(TRANSP.QS)));
        end
        adt=max(adt,1/365/24/60);    % use at least a 1 minute time step     
        if COAST.i_mc>1
            adt=min(adt,TIME.adt);   % make sure to use the minimum adt that is computed for all coastal segments
        end
        if TIME.it>1
            adt0=TIME.dt;
            adt=min(adt,10*adt0);     % make sure adt is not more than 10 times larger than previous cycle timestep 'adt0'
            adt=max(adt,0.01*adt0);   % make sure adt is not more than 100 times smaller than previous cycle timestep 'adt0'
        else
            adt=min(adt,1/365);       % use at most 1 day time step initially   
        end

    end
    TIME.adt=adt;
end


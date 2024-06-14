function [OP]=overwash_potential(T,Ho,hc,tanbeta)
%function [OP]=overwash_potential(T,Ho,hc,tanbeta)
%
% INPUT: 
%      T        : peak wave period        
%      Ho       : wave height
%      hc       : barrier crest elevation (now used as S.Bheight) 
%      tanbeta  : beach slope
% OUTPUT:
%      OP       : overwash potential
% function 
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

    OP=logical(zeros);
    g = 9.81;
    Lo =g*T.^2./2./pi; % the deep-water wave length
    % tanbeta : beach slope
    % hc : barrier crest elevation (now used as S.Bheight)
    z=0;  % water level
    OP=1.1*(0.35*tanbeta*sqrt(Ho)*Lo+sqrt(Ho*Lo*0.563*(tanbeta)^2+0.004)/2)+z-hc;
    % OP==1.1*(0.35*tanbeta*sqrt(Ho)*Lo+sqrt(Ho*Lo*0.563*(tanbeta)^2+0.004)/2)+z-hc>1;
    
end

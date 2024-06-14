function [prob]=get_dirspr_prob(dspr,dirout,mu)
% function [prob]=get_dirspr_prob(dspr,dirout,mu)
% 
% INPUT:
%     dspr    directional spreading [°]
%     dirout  (optional) vector with directions w.r.t. main direction (with default of [-40:10:40])
%     mu      (optional) main wave direction, used to correct the dirout (directions used = dirout-mu)
%
% OUTPUT:
%     prob    vector with probabilities with wave directions given the directional spectrum
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

    if nargin<3
        mu=0;
    end
    sigma=dspr; 
    if nargin<2
        ddir=10;
        dirout=[-40:10:40];
    end
    ddir=diff(dirout);
    
    prob=[];
    for pp=1:length(dirout)
        ddir1=ddir(max(pp-1,1));
        ddir2=ddir(min(pp,length(ddir)));
        p1=1/(sigma*sqrt(2*pi)).*exp((-(dirout(pp)-ddir1/2-mu).^2)/(2*sigma.^2));
        p2=1/(sigma*sqrt(2*pi)).*exp((-(dirout(pp)+ddir2/2-mu).^2)/(2*sigma.^2));    
        prob(pp)=(ddir1+ddir2)/2*(p2+p1)/2; % compute probability of area underneath normal distribution for 'meandir-ddir/2' to 'meandir+ddir/2'
    end
end

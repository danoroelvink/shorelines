function [kd]=wave_diffraction_coeff(omega,kdform,wdform,dspr)
% function [kd]=wave_diffraction_coeff(omega,kdform,wdform,dspr)
%
% The wave height reduction as a result of wave diffraction is computed in this routine. 
% 
% INPUT: 
%    omega   : angle difference in degrees [°]
%    kdform  : method for computing wave height attenuation due to diffraction (either 'Kamphuis' or 'Roelvink')
%    wdform  : method for computing impact of directional spreading on direction of diffracted waves (either 'Roelvink' or 'Dabees'/'Kamphuis')
%    dpsr    : (optional) directional spreading [°]
%
% OUTPUT:
%    kd      : wave height reduction factor [-]
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

    % create discretized components of the directional spectrum with their probabilities
    kd=[];
    for j=1:size(omega,1)
        
        % determine directional distribution probability per directional sector
        if strcmpi(wdform,'Roelvink')
            dirout=[0];
            omr=omega(j,:);
            prob=1;
        else
            dirout=[-50:5:50];
            omr=repmat(dirout(:),[1,length(omega(j,:))])+repmat(omega(j,:),[length(dirout),1]);
            prob=get_dirspr_prob(dspr,dirout);
            prob=prob/sum(prob);
        end
        
        % determine wave height reduction factor due to diffraction
        if kdform=='Roelvink'

            fac1=0.5;
            pwr=4;
            om_x=(omr+90)/180;
            kdi=1-exp(-abs((fac1./om_x).^pwr));
            kdi(omr<-90)=1;

        elseif kdform=='Kamphuis'

            omk=-omr;   % Note: we define omega positive towards shadow zone
            kdi=ones(size(omk));           
            %kdi(omk < -90 ) = 0;
            %kdi(omk <= 0 & omk >=-90) = 0.71+ 0.0093*omk(omk <= 0 & omk >=-90) + 0.0000156*omk(omk <= 0 & omk >=-90); % formulation that goeas to 0 at 90° and connects better at 0°
            kdi(omk < -86.25 ) = 0;
            kdi(omk <= 0 & omk >=-86.25) = max(0.69+0.008*omk(omk <= 0 & omk >=-86.25),0); 
            kdi(omk > 0 & omk <= 40) = 0.71+0.356686*sind(omk(omk >  0 & omk <= 40));
            kdi(omk > 40 & omk <= 90) = min(0.83+0.17*sind(omk(omk > 40 & omk <= 90)),1);
        end

        % recombine the discretized components of the directional spectrum based on their energy (quadratic)
        kdiprob=kdi.^2.*repmat(prob(:),[1,size(kdi,2)]);
        kd(j,:)=(sum(kdiprob,1)).^0.5;
    end
end            

function [H,Dw,Urms,k,C,Cg,theta,Fy,vw] = wave_cur_1D(H0,T,theta0,alpha,gamma,cf,x,h,n)
% function [H,Dw,Urms,k,C,Cg,theta,Fy,vw] = wave_cur_1D(H0,T,theta0,alpha,gamma,cf,x,h,n)
% 
% WAVE_CUR_1D Solve Snel's Law, wave energy balance and wave driven current
% on a 1D profile, for a range of water levels.
%
% INPUT:
%     H0        : Hrms wave height at seaward boundary (m)
%     T         : Peak wave period (s)
%     theta0    : offshore wave direction w.r.t.coast normal (deg)
%     alpha     : wave dissipation coefficient (-)
%     gamma     : wave breaking coefficient (-)
%     cf        : bed friction coefficient
%     x         : cross-shore distance (row vector) (m)
%     h         : water depth, 2D matrix [nx nT] (m)
%
% OUTPUT:
%     H         : Hrms wave height, 2D matrix [nx nT] (m)
%     Dw        : wave dissipation, 2D matrix [nx nT] (W/m2)
%     Urms      : rms orbital velocity, 2D matrix [nx nT] (m/s)
%     k         : wave number, 2D matrix [nx nT] (1/m)
%     C         : wave celerity, 2D matrix [nx nT] (m/s)
%     Cg        : wave group speed, 2D matrix [nx nT] (m/s)
%     theta     : wave angle w.r.t. coast normal, 2D matrix [nx nT] (deg)
%     Fy        : longshore component of wave force, 2D matrix [nx nT] (N/m2)
%     vw        : longshore wave-driven current velocity, 2D matrix [nx nT] (m/s)
% 
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2024 IHE Delft & Deltares
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

    g=9.81;
    rho=1025;
    eps=1e-6;
    if ~isempty(n)
    cf=9.81*(n^2)./h.^(1/3);
    end
    dx=x(2)-x(1);
    
    % Solve dispersion relation
    [k,C,Cg]=get_disper(h,T);
    Hmax=0.88./k.*tanh(gamma*k.*h/0.88);
    Emax=0.125*rho*g*Hmax.^2;
    
    % Solve wave direction from Snel's Law
    try
        theta=asind(sind(theta0).*C./C(1,:)); 
    catch
        theta=theta0;
    end

    E=zeros(size(h));
    beta=zeros(size(h));
    E(1,:)=max(1/8*rho*g.*H0.^2,eps);
    for i=1:length(x)-1
        beta(i,:)=2*alpha./T*dx.*exp(-Emax(i,:)./abs(E(i,:)));
        E(i+1,:)=(E(i,:).*Cg(i,:).*cosd(theta(i,:))-beta(i,:).*Emax(i,:))./(Cg(i+1,:).*cosd(theta(i+1,:))+beta(i,:));
        beta(i,:)=2*alpha./T*dx.*exp(-Emax(i,:)./abs(E(i,:)+E(i+1,:))*2);
        E(i+1,:)=(E(i,:).*Cg(i,:).*cosd(theta(i,:))-beta(i,:).*Emax(i,:))./(Cg(i+1,:).*cosd(theta(i+1,:))+beta(i,:));
        E(i+1,:)=max(min(E(i+1,:),0.125*rho*g*(gamma*h(i+1,:)).^2),eps);
    end

    H=sqrt(8*E/rho/g);
    Dw=beta/dx.*(Emax+E);
    Fy=Dw./C.*sind(theta);
    Urms=1/sqrt(2)*pi*min(H./T,1)./sinh(k.*h);
    aa=1./Urms.^2;bb=1.16^2;cc=-(Fy/rho./cf./Urms).^2;
    vw=sqrt((-bb+sqrt(bb.^2-4*aa.*cc))/2./aa).*sign(theta);
    Urms(isnan(Urms))=0;
    vw(isnan(vw))=0; 
    % vw=Fy/rho./cf./Urms;
    % vw=sqrt(abs(Fy)/rho./cf).*sign(Fy);
end

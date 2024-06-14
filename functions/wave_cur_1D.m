function [H,Dw,Urms,k,C,Cg,theta,Fy,vw] = wave_cur_1D(H0,T,theta0,alpha,gamma,Cf,x,h);
% WAVE_CUR_1D Solve Snel's Law, wave energy balance and wave driven current
%             on 1D profile, for a range of water levels.
%% Input
%  H0         : Hrms wave height at seaward boundary (m)
%  T          : Peak wave period (s)
%  theta0     : offshore wave direction w.r.t.coast normal (deg)
%  alpha      : wave dissipation coefficient (-)
%  gamma      : wave breaking coefficient (-)
%  Cf         : bed friction coefficient
%  x          : cross-shore distance (row vector) (m)
%  h          : water depth, 2D matrix [nx nT] (m)
%% Output
%  H          : Hrms wave height, 2D matrix [nx nT] (m)
%  Dw         : wave dissipation, 2D matrix [nx nT] (W/m2)
%  Urms       : rms orbital velocity, 2D matrix [nx nT] (m/s)
%  k          : wave number, 2D matrix [nx nT] (1/m)
%  C          : wave celerity, 2D matrix [nx nT] (m/s)
%  Cg         : wave group speed, 2D matrix [nx nT] (m/s)
%  theta      : wave angle w.r.t. coast normal, 2D matrix [nx nT] (deg)
%  Fy         : longshore component of wave force, 2D matrix [nx nT] (N/m2)
%  vw         : longshore wave-driven current velocity, 2D matrix [nx nT] (m/s)

g=9.81;
rho=1025;
dx=x(2)-x(1);
%tic
% Solve dispersion relation
[k,C,Cg]=get_disper(h,T);
Hmax=0.88./k.*tanh(gamma*k.*h/0.88);
Emax=0.125*rho*g*Hmax.^2;

% Solve wave direction from Snel's Law
try
theta=asind(sind(theta0).*C./C(1,:)); 
catch
    theta0
end

E=zeros(size(h));
beta=zeros(size(h));
E(1,:)=1/8*rho*g*H0^2;
for i=1:length(x)-1
   beta(i,:)=2*alpha/T*dx*exp(-Emax(i,:)./E(i,:));
   E(i+1,:)=(E(i,:).*Cg(i,:).*cosd(theta(i,:))-beta(i,:).*Emax(i,:))./(Cg(i+1,:).*cosd(theta(i+1,:))+beta(i,:));
   beta(i,:)=2*alpha/T*dx*exp(-Emax(i,:)./(E(i,:)+E(i+1,:))*2);
   E(i+1,:)=(E(i,:).*Cg(i,:).*cosd(theta(i,:))-beta(i,:).*Emax(i,:))./(Cg(i+1,:).*cosd(theta(i+1,:))+beta(i,:));
   E(i+1,:)=max(min(E(i+1,:),0.125*rho*g*(gamma*h(i+1,:)).^2),eps);
end

H=sqrt(8*E/rho/g);
Dw=beta/dx.*(Emax+E);
Fy=Dw./C.*sind(theta);
Urms=1/sqrt(2)*pi*H./T./sinh(k.*h);
aa=1./Urms.^2;bb=1.16^2;cc=-(Fy/rho./Cf./Urms).^2;
vw=sqrt((-bb+sqrt(bb.^2-4*aa.*cc))/2./aa).*sign(theta);
% vw=Fy/rho./Cf./Urms;
%vw=sqrt(abs(Fy)/rho./Cf).*sign(Fy);
%toc



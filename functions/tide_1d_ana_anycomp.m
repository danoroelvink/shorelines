function [h,vt,surfslope] = tide_1d_ana_anycomp(eta,detads, ...
                  phi,vmean,surfslope,Ttide,nT,k, ...
                  Cf,hmin,x,zb)
%TIDE_1D Analytical solution of tidal and mean velocity over profile 
%% Input:
%  eta        : vertical amplitude of tide components(m)
%  detads     : longshore gradient of eta (-)
%  phi        : phase of vertical tide components(deg)
%  surfslope  : mean longshore surface slope (-)
%  Ttide      : tidal period (s)
%  nT         : number of time points in tidal cycle (-)
%  k          : alongshore wave number of tidal wave components (rad/m)
%  Cf         : bed friction coefficient (-)
%  hmin       : minimum water depth (m)
%  x          : cross-shore distance (row vector) (m)
%  zb         : bed level (row vector) (m)
%% Output
%  h          : 2D matrix [nx nt] with water depth profiles
%  vt         : 2D matrix [nx nt] with intratidal velocity profiles

%disp('computing tidal velocity profiles...');
g=9.81;
deg2rad=pi/180;
ncomp=length(eta);
phi=phi*deg2rad;
omega(1)=2*pi/Ttide;
for i=2:ncomp
    omega(i)=i*omega(1);
end
dT=Ttide/nT;
t=[dT:dT:Ttide];
etat=zeros(size(t));
for i=1:ncomp
    etat=etat+eta(i)*cos(omega(i)*t-phi(i));
end
tt=repmat(t,length(x),1);
h=etat-zb';
h=max(h,hmin);
hmean=mean(h,2);
%% Correction for gradient of tidal amplitude
amp=sqrt(sum(k.^2.*eta.^2+detads.^2));

%% Solve labda and vamp analytically2
%tic
aa=(3*pi/8*Cf/omega(1)./hmean).^2;bb=1;cc=-(g*amp/omega(1))^2;
vamp=max(sqrt((-bb+sqrt(bb^2-4*aa*cc))/2./aa),eps);
labda=3*pi/8*Cf*vamp;
%% Compute the coefficients a and b of the analytical solution for v due to tide
for i=1:ncomp
    ph(:,:,i)=labda./hmean/omega(i);
    a(:,:,i)= g/omega(i).*(        1./(1+ph(:,:,i).^2)*k(i)*eta(i)-ph(:,:,i)./(1+ph(:,:,i).^2)*detads(i));
    b(:,:,i)=-g/omega(i).*(ph(:,:,i)./(1+ph(:,:,i).^2)*k(i)*eta(i)+        1./(1+ph(:,:,i).^2)*detads(i));
end
%% Solution of tidal velocity for all x and t
vt=zeros(size(h));
for i=1:ncomp
    vt=vt +a(:,:,i).*cos(omega(i)*tt-phi(:,i))+b(:,:,i).*sin(omega(i)*tt-phi(:,i));
end
%vt=vt-g*h*surfslope./labda;
if length(x)==1
    surfslope=-vmean*labda/g/hmean;
end
vt=vt-mean(vt,2)-g*hmean*surfslope./labda;  
%toc


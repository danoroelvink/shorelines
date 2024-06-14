function [Slong]=soulsby_van_rijn(h,T,k,H,v,ks,hmin,D50,D90)
%
% Compute sand transport according to Soulsby - van Rijn
% ------------------------------------------------------
g=9.81;delta=1.65;rho=1025;nu=1.e-6;kappa=0.39;Acal=0.2;
z0=ks/30;
%z0=0.006;
Dstar=(g*delta/nu^2)^(1/3)*D50;
dry=h<hmin;
h=max(h,hmin);
Urms=1/sqrt(2)*pi*H./T./sinh(k.*h);
if D50<=0.5e-3
    Ucr=0.19*D50^0.1*log10(4*h/D90);
else
    Ucr=8.5*D50^0.6*log10(4*h/D90);
end
Cf=(kappa./(log(h/z0)-1)).^2;
umod=sqrt(v.^2+0.018./Cf.*Urms.^2);
ksi=(umod-Ucr).^2.4;
ksi(umod<Ucr)=0;
Asb=0.005*h.*(D50./h/delta/g/D50).^1.2;
Ass=0.012*D50*Dstar^(-0.6)/(delta*g*D50)^1.2;
Sby=Acal*Asb.*v.*ksi;
Ssy=Acal*Ass.*v.*ksi;
Stoty=Sby+Ssy;
%% Set values to 0 at dry points
Urms(dry)=0;
Ucr(dry)=0;
Cd(dry)=0;
ksi(dry)=0;
Asb(dry)=0;
Ass(dry)=0;
Sbx(dry)=0;
Sby(dry)=0;
Ssx(dry)=0;
Ssy(dry)=0;
Stotx(dry)=0;
Stoty(dry)=0;
Slong=Stoty*3600*24*365;

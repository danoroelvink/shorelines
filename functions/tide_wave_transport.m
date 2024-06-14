function [Slongtot,Slong_mean,Slong,v,vt,vw,Hrms,h] = tide_wave_transport ...
         (eta,detads,phi,surfslope,Ttide,nT,ktide,Cf,hmin, ...
          Hrms0,Tp,theta0,alpha,gamma, ...
          ks,D50,D90,hclosure, ...
          x,zb);
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%% Compute tidal velocity, intra-tide and residual
[h,vt] = tide_1d_ana_anycomp(eta,detads,phi,0,surfslope,Ttide,nT,ktide,Cf,hmin,x,zb);
%% Compute wave decay throughout the tide and wave-driven current 
[Hrms,Dw,Urms,k,C,Cg,theta,Fy,vw] = wave_cur_1D(Hrms0,Tp,theta0,alpha,gamma,Cf,x,h);
%% Total current (x,t)
v=vt+vw;
%% Longshore transport for each x and t
[Slong]=soulsby_van_rijn(h,Tp,k,Hrms,v,ks,hmin,D50,D90);
%% Average longshore transport profile
Slong(h==hmin)=0;
Slong(h>hclosure)=0;
Slong_mean=mean(Slong,2);
hmean=mean(h,2);
%% Total longshore transport
Slongtot=sum(Slong_mean)*(x(2)-x(1));
end


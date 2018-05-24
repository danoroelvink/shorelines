function [kh,c] = GUO2002(Tp,h)
% function kh=GUO2002(Tp,h) computes the kh number of the waves for shallow water
%
% INPUT:
%     Tp       peak wave period [s]
%     h        water depth [m]
%
% OUTPUT:
%     kh       wave number x water depth
%
% by: B.J.A.Huisman, 2017 (Deltares)

omega        = (2*pi)./Tp;
kh           = omega.^2.*h/9.81.*(1-exp(-1*(omega.*(h/9.81).^0.5).^2.5)).^-0.4;
c            = omega.*(h./kh);
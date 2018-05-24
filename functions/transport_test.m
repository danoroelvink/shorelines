clc 
clear all
%function [QS,adt,philoc]=transport2(phic,philoc,phiw,S,it,s)
%% Transform to nearshore
S.d=15;
S.Hso=1.4;
S.ds0=250;
phiw=[-90:0.01:90]*pi/180;
phic=0;
philoc=atan2(sin(phic-phiw),cos(phic-phiw));
S.tper=6;
S.ddeep=25;
S.dnearshore=8;
S.phif=0*pi/180;                                                         % Orientation of the foreshore [°N] <- only relevant for 'KAMP', 'MILH' or 'VR14'
S.trform='CERC';                                                           % switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
S.b=0.5e6;                                                                   % CERC : coeff in simple cerc formula
S.b2=1.28e6;
S.tper=6;                                                                  % KAMP & MILH : peak wave period [s]
S.d50=2.0e-4;                                                              % KAMP & MILH & VR14 : median grain diameter [m]
S.porosity=0.4;                                                            % KAMP & MILH & VR14 : S.porosity (typically 0.4) [-]
S.tanbeta=0.03;                                                            % KAMP & MILH & VR14 : mean bed slope [ratio 1/slope]
S.rhos=2650;                                                               % KAMP & MILH & VR14 : density of sand [kg/m3]
S.rhow=1025;                                                               % KAMP & MILH & VR14 : density of water [kg/m3]
S.g=9.81;                                                                  % KAMP & MILH & VR14 : gravitational acceleration [m2/s]
S.alpha=1.8;                                                               % KAMP & MILH & VR14 : calibration factor for point of breaking (S.alpha = 1.8 for Egmond data)
S.gamma=0.72;                                                              % KAMP & MILH & VR14 : breaking coefficient (Hs/h) with 5% breaking waves
S.Pswell=20; 


%if strcmpi(S.trform,'CERC2')
    %           phi9=[0:.01:70]*pi/180;id=find(cos(phi9).^(6/5).*sin(phi9)==max(cos(phi9).^(6/5).*sin(phi9)));
    %           thetacrit=thetacrit0*(phi9(id)*180/pi/45);
    [kh0,c0]=GUO2002(S.tper,S.ddeep);
    [kh1,c1]=GUO2002(S.tper,S.dnearshore);

    
%elseif ~strcmpi(S.trform,'CERC') && ~strcmpi(S.trform,'CERC2')
%     refraction on foreshore (S.ddeep -> S.dnearshore)
    if isempty(S.phif)
        phif=phic;
    end
    philoc0=atan2(sin(S.phif-phiw),cos(S.phif-phiw));
    [kh0,c0]=GUO2002(S.tper,S.ddeep);
    [kh1,c1]=GUO2002(S.tper,S.dnearshore);
    philoc1=asin((c1./c0).*sin(philoc0));
    
    % refraction in the nearshore
    philoc1B = philoc1+(phic-S.phif);
    hbr=real(((S.Hso.^2 .* c1 .* cos(philoc))./(S.alpha .* S.gamma.^2 .* S.g.^0.5)).^0.4);
    hsbr=hbr*S.gamma;
    [kbr,cbr]=GUO2002(S.tper,hbr);
    sphibr=asin((cbr./c1).*sin(philoc1B));
    
[kh0,c0]=GUO2002(S.tper,S.ddeep);
    [kh1,c1]=GUO2002(S.tper,S.dnearshore);
    hbr=real(((S.Hso.^2 .* c1 .* cos(philoc))./(S.alpha .* S.gamma.^2 .* S.g.^0.5)).^0.4);
    hsbr=hbr*S.gamma;
    [kbr,cbr]=GUO2002(S.tper,hbr);
    phib=asin((cbr./c0).*sin(phiw));
    philocB=atan2(sin(phic-phib),cos(phic-phib));
    phib=-phib;
% end



% if it==0   %four boundary conditions
%     QS=zeros(size(s));
%     S.QS_start=QS(1);
%     S.QS_end=QS(end);
% end










%% Transport : CERC
%if strcmpi(S.trform,'CERC')
    k=0.2;                                                      % using CERC (1984) value of k(SPM,Hs)=0.39
    % b_theory = k * (S.rhow * S.g^0.5 / (16 * sqrt(S.gamma)* (S.rhos-S.rhow) * (1-S.porosity)));  % b_theory = 0.0946 * 365*24*60*60 = 2.9833E+6
    QS=S.b*S.Hso^2.5*sin(2*(philoc));QScerc=QS;
    %Smax=S.b*S.Hso^2.5;
    adt=0.25*S.ds0^2/(S.b*S.Hso^2.5*1)*S.d %(365*24*60*60)^2 %assumption cos2deltb=1
    
%end
%if strcmpi(S.trform,'CERC2')
    k=0.12;                                                      % using CERC (1984) value of k(SPM,Hs)=0.39
    b1 = k * (S.rhow * S.g^0.5 / (16 * sqrt(S.gamma)* (S.rhos-S.rhow) * (1-S.porosity)));  % b_theory = 0.0946 * 365*24*60*60 = 2.9833E+6
    b2 = b1 * ((S.gamma.*S.g).^0.5 /(2*pi)).^0.2;  % b_theory = 0.0946 * 365*24*60*60 = 2.9833E+6
    QS=365*24*60*60*b2*S.Hso.^(12/5)*S.tper.^(1/5).*cos(philoc).^(6/5).*sin(philoc);QScerc2=QS;
    %Smax=S.b*S.Hso^2.5;
    adt=(0.5*(S.ds0.^2))./((b2/S.d*S.Hso.^(12/5)).*S.tper.^(1/5))/(365*24*60*60);%*cos(philoc_cor).^(1/5).*((cos(philoc_cor).^2)-(6/5).*(sin(philoc_cor).^2)))/(365*24*60*60);
    %adt=min(adt(adt > 0))
    % adt=real(adt)
    
%end

%% Transport : Kamphuis
%if strcmpi(S.trform,'KAMP')
QSkampmass=2.33 * S.rhos/(S.rhos-S.rhow) .* S.tper.^1.5 .* S.tanbeta.^0.75 .* S.d50.^-0.25 .* hsbr.^2 .* (abs(sin(2*sphibr)).^0.6.*sign(sphibr) );    
%QSkampmass=2.33 * S.rhos/(S.rhos-S.rhow) .* S.tper.^1.5 .* S.tanbeta.^0.75 .* S.d50.^-0.25 .* hsbr.^2 .* (abs(sin(2*phib)).^0.6.*sign(phib) );
    QS = 365*24*60*60*(QSkampmass /(S.rhos-S.rhow)) /(1.0-S.porosity);QSkamp=QS;
    adt=(0.5*(S.ds0.^2))./(QS).*S.d.*phib;
    adt=min(adt(adt > 0))
%     %% Transport : Mil-Homens
% %elseif strcmpi(S.trform,'MILH')
%     QSmilhmass=0.15*S.rhos/(S.rhos-S.rhow) .* S.tper.^0.89 .* S.tanbeta.^0.86 .* S.d50.^-0.69 .* hsbr.^2.75 .* real((sin(2*sphibr)).^0.5);
%     QS = 365*24*60*60*(QSmilhmass /(S.rhos-S.rhow)) /(1.0-S.porosity);QSmilh=QS;
%     %% Transport : Van Rijn (2014)
% %elseif strcmpi(S.trform,'VR14')
%     kswell=0.015*S.Pswell+(1-0.01*S.Pswell);
%     vwave=0.3*real(sin(2*sphibr)).*(S.g.*hsbr).^0.5;
%     vtide=0;
%     vtotal=vwave+vtide;
%     QSvr14mass=0.0006 .* kswell .* S.rhos .* S.tanbeta.^0.4 .*S.d50.^-0.6 .* hsbr.^2.6 .* abs(vtotal);
%     QS = 365*24*60*60*(QSvr14mass /(S.rhos-S.rhow)) /(1.0-S.porosity);QSvr14=QS;
    %%Transport: CERC (vitousek,Barnard2015)
%elseif strcmpi(S.trform,'CERC3')
    QS=S.b2*hsbr.^2.5.*sin(2*(philocB));QScerc3=QS;
    adt=0.25*S.ds0.^2./(S.b*hsbr.^2.5)*S.d;
    adt=min(adt)
%end
% 
phiw=phiw*180/pi;
% plot(phiw,QScerc,phiw,QScerc2,phiw,QSkamp,phiw,QScerc3,'LineWidth',2)
% legend('CERC','CERC2','KAMP','CERC3')
% % 
% % Qm=max(QScerc);
% % thtcr=interp1(QScerc,phiw,Qm) 
% 
%  ssss=max(QScerc)
  QST={QScerc,QScerc2,QSkamp,QScerc3}
% for i = 1:4
% [QSmax(i),ind]=max(QST{i})
% thetacrit(i)=phiw(ind)
% end
% 
for i = 3
[QSmax(i),ind]=max(QST{i})
thetacritb(i)=sphibr(ind)*180/pi
end
% 
for i = 4
[QSmax(i),ind]=max(QST{i})
thetacritb(i)=philocB(ind)*180/pi
end

%plot (sphibr*180/pi,QSkamp,'--')

% plot(sphibr*180/pi,cos(2*sphibr));
plot (philocB*180/pi,QScerc3,'--')
plot (sphibr*180/pi,QSkamp,philocB*180/pi,QScerc3,'LineWidth',2)
legend('KAMP','CERC3')
%  print('QSequaton2','-dpng', '-r300')
 xlabel('Local wave angle (Breaking) [deg]');
 ylabel('Longshore transport [m^3/yr]');
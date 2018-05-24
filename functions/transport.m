function [QS,QSmax,adt]=transport(S,it,s,philoc_cor,hsbr,sphibr)
%% Transport : CERC
if strcmpi(S.trform,'CERC')
    k=0.2;                                                      % using CERC (1984) value of k(SPM,Hs)=0.39
    % b_theory = k * (S.rhow * S.g^0.5 / (16 * sqrt(S.gamma)* (S.rhos-S.rhow) * (1-S.porosity)));  % b_theory = 0.0946 * 365*24*60*60 = 2.9833E+6
    QS=S.b*S.Hso.^2.5*sin(2*(philoc_cor));QScerc=QS;
%         QS=(31536000)*sin(2*(philoc_cor));QScerc=QS; % for test1
    %Smax=S.b*S.Hso^2.5;
%     adt=0.25*S.ds0^2/(S.b*S.Hso^2.5*1)*S.d %(365*24*60*60)^2 %assumption cos2deltb=1
    QSmax=S.b*max(S.Hso)^2.5*sin(2*(45*pi/180));
%     QSmax=(31536000)*sin(2*(45*pi/180)); %for test1
end
if strcmpi(S.trform,'CERC2')
    k=0.1;                                                      % using CERC (1984) value of k(SPM,Hs)=0.39
    b1 = k * (S.rhow * S.g^0.5 / (16 * sqrt(S.gamma)* (S.rhos-S.rhow) * (1-S.porosity)));  % b_theory = 0.0946 * 365*24*60*60 = 2.9833E+6
    b2 = b1 * ((S.gamma.*S.g).^0.5 /(2*pi)).^0.2;  % b_theory = 0.0946 * 365*24*60*60 = 2.9833E+6
    QS=365*24*60*60*b2*S.Hso.^(12/5)*S.tper.^(1/5).*cos(philoc_cor).^(6/5).*sin(philoc_cor);QScerc2=QS;
    %Smax=S.b*S.Hso^2.5;
%     adt=(0.5*(S.ds0.^2))./((b2/S.d*S.Hso.^(12/5)).*S.tper.^(1/5))/(365*24*60*60)%*cos(philoc_cor).^(1/5).*((cos(philoc_cor).^2)-(6/5).*(sin(philoc_cor).^2)))/(365*24*60*60);
    %adt=min(adt(adt > 0))
    % adt=real(adt)
    thetacrit=42.39*pi/180;
    QSmax=365*24*60*60*b2*max(S.Hso).^(12/5)*S.tper.^(1/5).*cos(42.39*pi/180).^(6/5).*sin(42.39*pi/180);
end

%% Transport : Kamphuis
if strcmpi(S.trform,'KAMP')
    QSkampmass=2.33 * S.rhos/(S.rhos-S.rhow) .* S.tper.^1.5 .* S.tanbeta.^0.75 .* S.d50.^-0.25 .* hsbr.^2 .* (abs(sin(2*sphibr)).^0.6.*sign(sphibr));
    QS = 365*24*60*60*(QSkampmass /(S.rhos-S.rhow)) /(1.0-S.porosity);QSkamp=QS;
    phibmax=15.5289*pi/180;
    QSkampmassmax=2.33 * S.rhos/(S.rhos-S.rhow) .* S.tper.^1.5 .* S.tanbeta.^0.75 .* S.d50.^-0.25 .* hsbr.^2 .* (abs(sin(2*phibmax)).^0.6.*sign(phibmax));
    QSmax = 365*24*60*60*(QSkampmassmax /(S.rhos-S.rhow)) /(1.0-S.porosity);
    QSmax=max(QSmax);
%     adt=(0.5*(S.ds0.^2))./(QSmax).*S.d.*sphibr;
%     adt=min(adt(adt > 0))
    
    %% Transport : Mil-Homens
elseif strcmpi(S.trform,'MILH')
    QSmilhmass=0.15*S.rhos/(S.rhos-S.rhow) .* S.tper.^0.89 .* S.tanbeta.^0.86 .* S.d50.^-0.69 .* hsbr.^2.75 .* real((sin(2*sphibr)).^0.5);
    QS = 365*24*60*60*(QSmilhmass /(S.rhos-S.rhow)) /(1.0-S.porosity);QSmilh=QS;
    %% Transport : Van Rijn (2014)
elseif strcmpi(S.trform,'VR14')
    kswell=0.015*S.Pswell+(1-0.01*S.Pswell);
    vwave=0.3*real(sin(2*sphibr)).*(S.g.*hsbr).^0.5;
    vtide=0;
    vtotal=vwave+vtide;
    QSvr14mass=0.0006 .* kswell .* S.rhos .* S.tanbeta.^0.4 .*S.d50.^-0.6 .* hsbr.^2.6 .* abs(vtotal);
    QS = 365*24*60*60*(QSvr14mass /(S.rhos-S.rhow)) /(1.0-S.porosity);QSvr14=QS;
    %% Transport: CERC-breaking conditions (vitousek,Barnard2015)
elseif strcmpi(S.trform,'CERC3')
    QS=S.b*hsbr.^2.5.*sin(2*(sphibr));QScerc3=QS;
    
%     adt=0.25*S.ds0.^2./(S.b*hsbr.^2.5)*S.d;
%     adt=min(adt)
    philocBmax=16.5436*pi/180;
    QSmax=S.b*hsbr.^2.5.*sin(2*(philocBmax));
    QSmax=max(QSmax);
end


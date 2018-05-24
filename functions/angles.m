function [philoc,thetacrit,hsbr,sphibr,hbr]=angles(S,s,x,y,n,phiwi)
sphi=[];
phic=[];
for i=1:n
    sphi(i)=.5*(s(i)+s(i+1));
    phic(i)=2*pi-atan2(y(i+1)-y(i),x(i+1)-x(i));
end


% if S.wave_interaction
%     philoc=atan2(sin(phic-phiw_cd),cos(phic-phiw_cd));
% else
%     philoc=atan2(sin(phic-phiw),cos(phic-phiw));
% end
philoc=atan2(sin(phic-phiwi),cos(phic-phiwi));
% Define critical angles
if strcmpi(S.trform,'CERC')
    theta=pi/4;
elseif strcmpi(S.trform,'CERC2')
    theta=42.39*pi/180;
elseif strcmpi(S.trform,'KAMP')
    theta=38.05*pi/180;
elseif strcmpi(S.trform,'CERC3')
    theta=41.42*pi/180;
end


thetacrit0 = repmat(theta,[1 n]);
thetacrit=thetacrit0;

%% Transform to nearshore

if ~strcmpi(S.trform,'CERC') && ~strcmpi(S.trform,'CERC2')
    %refraction on foreshore (S.ddeep -> S.dnearshore)
    if isempty(S.phif)
        S.phif=phic;
    end
    philoc0=atan2(sin(S.phif-phiwi),cos(S.phif-phiwi));
    [kh0,c0]=GUO2002(S.tper,S.ddeep);
    [kh1,c1]=GUO2002(S.tper,S.dnearshore);
    philoc1=asin((c1./c0).*sin(philoc0));
    
    % refraction in the nearshore
    philoc1B = philoc1+(phic-S.phif);
    hbr=real(((S.Hso.^2 .* c1 .* cos(philoc))./(S.alpha .* S.gamma.^2 .* S.g.^0.5)).^0.4);
    hsbr=hbr*S.gamma;
    [kbr,cbr]=GUO2002(S.tper,hbr);
    sphibr=asin((cbr./c1).*sin(philoc1B));
else
    hsbr=[];
    sphibr=[];
    hbr=[];
end
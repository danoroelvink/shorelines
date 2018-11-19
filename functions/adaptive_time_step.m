% function [adt]=adaptive_time_step(x_mc,y_mc,S,phiw,Hg,phiwg)
function [adt]=adaptive_time_step(x_mc,y_mc,S,phiw)

n_mc=length(find(isnan(x_mc)))+1;
for i_mc=1:n_mc
    n_mc=length(find(isnan(x_mc)))+1; %test
    if i_mc>n_mc
        break
    end
    
    [s,x,y,x_mc,y_mc,ar]=make_sgrid_mc(x_mc,y_mc,S.ds0,i_mc,S.smoothfac);
    n=length(x)-1;
    
   dsmin=min(diff(s));   %% adaptive time step calculations  
    %% Transport : CERC
    if strcmpi(S.trform,'CERC')
        adtc(i_mc)=0.25*dsmin^2/(S.b*max(S.Hso)^2.5*1)*S.d; %(365*24*60*60)^2 %assumption cos2deltb=1
%         adtc(i_mc)=0.25*dsmin^2/(31536000)*S.d %fot test1
    end
    if strcmpi(S.trform,'CERC2')
        k=0.12;                                                      % using CERC (1984) value of k(SPM,Hs)=0.39
        b1 = k * (S.rhow * S.g^0.5 / (16 * sqrt(S.gamma)* (S.rhos-S.rhow) * (1-S.porosity)));  % b_theory = 0.0946 * 365*24*60*60 = 2.9833E+6
        b2 = b1 * ((S.gamma.*S.g).^0.5 /(2*pi)).^0.2;  % b_theory = 0.0946 * 365*24*60*60 = 2.9833E+6
        adtc(i_mc)=(0.5*(dsmin^2))./((b2/S.d*S.Hso.^(12/5)).*S.tper.^(1/5))/(365*24*60*60);%*cos(philoc_cor).^(1/5).*((cos(philoc_cor).^2)-(6/5).*(sin(philoc_cor).^2)))/(365*24*60*60);

        %adt=min(adt(adt > 0))
        % adt=real(adt)
    end
    
    %% Transport : Kamphuis
    if strcmpi(S.trform,'KAMP')
        QSkampmassmax=2.33 * S.rhos/(S.rhos-S.rhow) .* S.tper.^1.5 .* S.tanbeta.^0.75 .* S.d50.^-0.25 .* hsbr.^2 .* (abs(sin(2*15.5289*pi/180)).^0.6.*sign(15.5289*pi/180));
        QSmax = 365*24*60*60*(QSkampmassmax /(S.rhos-S.rhow)) /(1.0-S.porosity);
        QSmax=max(QSmax);
        sphibr=abs(min(sphibr));
        adtc(i_mc)=(0.5*(dsmin^2))./(QSmax).*S.d.*(15.5289*pi/180);
% adtc(i_mc)=(0.01*(dsmin^2)).*S.d;
        %adt=min(adt(adt > 0))
        
        %% Transport : Mil-Homens
    elseif strcmpi(S.trform,'MILH')
        
        
        %% Transport : Van Rijn (2014)
    elseif strcmpi(S.trform,'VR14')
        
        
        %% Transport: CERC-breaking conditions (vitousek,Barnard2015)
    elseif strcmpi(S.trform,'CERC3')
        %hsbr=min(hsbr);
        adtcc=0.25*dsmin^2./(S.b*S.Hso.^2.5)*S.d;
%         adtcc=0.25*dsmin^2/(31536000)*S.d; %fot test1
        adtc(i_mc)=min(adtcc);
    end
    
end

adt=min(adtc);

function [phiw,S]=introduce_wave(S,WVC,WC,tnow)
    if ~isempty(S.WVCfile)
        %WVC.timenumi=timenuma;  %follow the consequences
        WVC.Hsi=interp1(WVC.timenum,WVC.Hs,tnow); %If the interpolation needs more computational .. We could interpolate in between
        WVC.Tpi=interp1(WVC.timenum,WVC.Tp,tnow);
        WVC.Diri=mod(atan2(interp1(WVC.timenum,sin(WVC.Dir*pi/180),tnow),interp1(WVC.timenum,cos(WVC.Dir*pi/180),tnow))*180/pi,360);
        S.Hso=interpNANs(WVC.Hsi);
        S.tper=interpNANs(WVC.Tpi);
        phiw=interpNANsDIR(WVC.Diri)*pi/180;
        %figure;plot(WVC.timenum,WVC.Dir,'b.',WVC.timenumi,WVC.Diri,'k.');
        % pause
        
    elseif ~isempty(S.Waveclimfile)
        iwc=round((rand*length(WC.Hs)+0.5));
        phiw=WC.dir(iwc)*pi/180;
        S.Hso=WC.Hs(iwc);
        S.tper=WC.Tp(iwc);
    else
        phiw=S.phiw0+(rand-.5)*pi/180*S.spread;
        S.Hso=S.Hso_in;
    end
%     if S.wave_interaction
%         [ Hg, phiwg ] = get_interpolated_wavefield( S.xg,S.yg,S.Hg_all,S.phiwg_all,S.Hso,phiw*180/pi);
%     end
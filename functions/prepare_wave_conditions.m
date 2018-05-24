function [WVC,WC,timenum0]=prepare_wave_conditions(S)
if ~isfield(S,'timenum0')
    if ~isempty(S.reftime)
        timenum0 = datenum(S.reftime,'yyyy-mm-dd'); %HH:MM:SS
    else
        timenum0 = 0;
    end
else
    timenum0=S.timenum0;
end
S.times(1)=timenum0;
    WVC=struct;
WC=struct;
if ~isempty(S.WVCfile)
    if S.RWSfiletype
        WVCraw=load(S.WVCfile);warning off;
%         WVC=struct;
        WVC.timenum=datenum([num2str(WVCraw(:,1))],'yyyymmddHHMM'); %,'%08.0f'),num2str(WVCraw(:,2),'%06.0f')
        WVC.Hs=interpNANs(WVCraw(:,2));
        WVC.Tp=interpNANs(WVCraw(:,3));
        WVC.Dir=interpNANsDIR(WVCraw(:,4));
    else
        WVCraw=load(S.WVCfile);warning off;
        WVC=struct;
        WVC.timenum=datenum([num2str(WVCraw(:,1),'%08.0f'),num2str(WVCraw(:,2),'%06.0f')],'yyyymmddHHMMSS');
        WVC.Hs=interpNANs(WVCraw(:,2));
        WVC.Tp=interpNANs(WVCraw(:,3));
        WVC.Dir=interpNANsDIR(WVCraw(:,5));
    end
   % timenuma=WVC.timenum;
  
elseif ~isempty(S.Waveclimfile)
    WCraw=load(S.Waveclimfile);
    WC=struct;
    WC.Hs=WCraw(:,1);
    WC.Tp=WCraw(:,2);
    WC.dir=WCraw(:,3)+S.Wavecorr;
  %  timenuma=timenum0;
else
%    timenuma=timenum0;
end
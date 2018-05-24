function [QS]=Sources_Sinks(QS,S,tnow,x_mc,y_mc)
if  strcmpi(S.sources_sinks,'FUNC') 
  SSraw=load(S.SSfile);warning off;
   
  for is=1:length(SSraw(:,1))/2
    SS=struct;
    SS.x_ss=SSraw((2*is-1),1);
    SS.y_ss=SSraw(2*is,1);
    SS.timenum=datenum([num2str(SSraw(:,2))],'yyyy');
    SS.QS_SS=interpNANs(SSraw(:,is+2));
    S.QS_SRC_SNK=interp1(SS.timenum,SS.QS_SS,tnow);
    SS=find(abs(x_mc-SS.x_ss) < 0.55*S.ds0); %& abs(y_mc-SS.y_ss)<0.55*S.ds0 )
    QS(SS)=S.QS_SRC_SNK;
  end
end


if strcmpi(S.sources_sinks,'CONST') 
     SSraw=load(S.SSfile);warning off;
     for is=length(SSraw(:,1))
         SS=struct;
         SS.x_ss=SSraw(is,1);
         SS.y_ss=SSraw(is,2);
         S.QS_SRC_SNK=SSraw(is,3);
         SS=find(abs(x_mc-SS.x_ss) < 0.55*S.ds0); %& abs(y_mc-SS.y_ss)<0.55*S.ds0 )
         QS(SS)=S.QS_SRC_SNK;
     end
end
end
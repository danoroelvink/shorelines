function [QS]=upwind_correction(philoc,n,cyclic,S,thetacrit,shadowS,shadowS_h,QSmax,QS)

for i=1:n
    if cyclic
        im1=mod2(i-1,n);
        im2=mod2(i-2,n);
        ip1=mod2(i+1,n);
        ip2=mod2(i+2,n);
    else
        im1=max(i-1,1);
        im2=max(i-2,1);
        ip1=min(i+1,n);
        ip2=min(i+2,n);
    end
    %if abs(philoc(i))>pi/4&&philoc(im1)<pi/4&&philoc(im1)>0
    if philoc(i)>thetacrit(i)&&philoc(im1)<thetacrit(i)&& ...
            philoc(im1)>0&&~shadowS(im1)&& ...
            (isempty(shadowS_h)||(~isempty(shadowS_h)&&~shadowS_h(ip1)&&~shadowS_h(ip2)))
        QS(i)=QSmax;
        if S.twopoints
            QS(ip1)=0.5*QSmax;
            QS(ip2)=0;
        else
            %QS(ip1)=0;
        end
        %elseif abs(philoc(i))>pi/4&&philoc(ip1)>-pi/4&&philoc(ip1)<0
    elseif philoc(i)<-thetacrit(i)&&philoc(ip1)>-thetacrit(i)&& ...
            philoc(ip1)<0&&~shadowS(ip1)&& ...
            (isempty(shadowS_h)||(~isempty(shadowS_h)&&~shadowS_h(im1)&&~shadowS_h(im2)))
        QS(i)=-QSmax;
        if  S.twopoints
            QS(im1)=-0.5*QSmax;
            QS(im2)=0;
        else
            %QS(im1)=0;
        end
    end
end
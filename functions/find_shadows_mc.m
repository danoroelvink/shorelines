function [ xS,yS,shadowS ] = find_shadows( x,y,x_mc,y_mc,phiw,hard )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
n=length(x)-1;
if hard==1
    crit=0;
else
    crit=0;
end
if n==0
    xS=[];
    yS=[];
    shadowS=[];
else
    f=180/pi;
    xS=.5*(x(1:n)+x(2:n+1));
    yS=.5*(y(1:n)+y(2:n+1));
    len=5*hypot(max(x_mc)-min(x_mc),max(y_mc)-min(y_mc));
    
    for i=1:n
        xw=[xS(i)+1*sin(phiw),xS(i)+len*sin(phiw)];
        yw=[yS(i)+1*cos(phiw),yS(i)+len*cos(phiw)];
        P1=InterX([x_mc;y_mc],[xw;yw]);
        
        shadowS(i)=size(P1,2)>crit;
        if 0
            figure(10)
            plot(x,y,x_mc,y_mc,xw,yw,'.-',P1(1,:),P1(2,:),'ro','linewidth',2);
            axis equal
            drawnow
            pause
        end
    end
end

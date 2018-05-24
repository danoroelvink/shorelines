function [ x,y, spit,width ] = find_spit_width_v2( s,x,y,phiw,shadow,spit_width )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here4
crit_width=spit_width;
eps=1;
width=zeros(size(x));
spit=logical(zeros(size(x)));
n=length(x)-1;
f=180/pi;
len=hypot(max(x)-min(x),max(y)-min(y));

for i=2:n-1
    % construct line normal to coast
    %     xl=[x(i)-eps*(y(i+1)-y(i-1))/hypot(x(i+1)-x(i-1),y(i+1)-y(i-1)),x(i)+len*(y(i+1)-y(i-1))/hypot(x(i+1)-x(i-1),y(i+1)-y(i-1))];
    %     yl=[y(i)+eps*(x(i+1)-x(i-1))/hypot(x(i+1)-x(i-1),y(i+1)-y(i-1)),y(i)-len*(x(i+1)-x(i-1))/hypot(x(i+1)-x(i-1),y(i+1)-y(i-1))];
    xl=[x(i)-eps*sin(phiw),x(i)+len*sin(phiw)];
    yl=[y(i)-eps*cos(phiw),y(i)+len*cos(phiw)];
    if i>1
        
        P1=InterX([x;y],[xl;yl]);
        P1=unique(P1','rows','stable')';
        spit(i)=size(P1,2)==2;
        if size(P1,2)==2
            width(i)=hypot(P1(1,1)-P1(1,2),P1(2,1)-P1(2,2));
        end
        if 0
            figure(10)
            if size(P1,2)==2
                plot(x,y,xl,yl,'ro-',P1(1,1),P1(2,1),'*k',P1(1,2),P1(2,2),'*b');
                title(['width=',num2str(width(i))]);
                P1
                pause
            end
        end
    else
        spit(i)=logical(0);
    end
    %     if i<n
    %         xx2=[x(i:n+1)];
    %         yy2=[y(i:n+1)];
    %         P2=InterX([xx2;yy2],[xl;yl]);
    %         spit(i)=spit(i)|size(P2,2)>1;
    %         if size(P2,2)>2
    %             width(i)=hypot(x(i)-P2(1,2),y(i)-P2(2,2));
    %         end
    %     end
    if width(i)>.1 & width(i)<crit_width & size(P1,2)==2 & shadow(i)==0
        if hypot(x(i)-P1(1,1),y(i)-P1(2,1))>hypot(x(i)-P1(1,2),y(i)-P1(2,2))
            [mindist, ip] = min(hypot(x-P1(1,1),y-P1(2,1)));
        else
            [mindist, ip] = min(hypot(x-P1(1,2),y-P1(2,2)));
        end
        
        if ip>1&&ip<n
            yesplot=0;
            if yesplot
                figure(10)
                %    plot(x,y,xx2,yy2,xl,yl,'.-',P1(1,:),P1(2,:),'ro',P2(1,:),P2(2,:),'bo','linewidth',2);
                plot(x,y,xl,yl,'.-',P1(1,:),P1(2,:),'ro',P1(1,2),P1(2,2),'k+','linewidth',2);
                axis equal
                %ylim([8000 12000]);
                title(['width ',num2str(round(width(i)))])
                hold on
                plot(x(ip),y(ip),'ok',x(ip-1),y(ip-1),'ob',x(ip+1),y(ip+1),'ob')
            end
            %% shift landward side of spit
            dn=min(crit_width-width(i),.1*(s(2)-s(1)));
            dx=-dn*(y(i+1)-y(i-1))/hypot(x(i+1)-x(i-1),y(i+1)-y(i-1));
            dy= dn*(x(i+1)-x(i-1))/hypot(x(i+1)-x(i-1),y(i+1)-y(i-1));
            x(i)=x(i)+dx;
            y(i)=y(i)+dy;
            %% shift seaward side of spit
            dn=-dn;
            dx=-dn*(y(ip+1)-y(ip-1))/hypot(x(ip+1)-x(ip-1),y(ip+1)-y(ip-1));
            dy= dn*(x(ip+1)-x(ip-1))/hypot(x(ip+1)-x(ip-1),y(ip+1)-y(ip-1));
            x(ip)=x(ip)+dx;
            y(ip)=y(ip)+dy;
            if yesplot
                plot(x,y,'--k')
                hold off
                drawnow
                pause
            end
        end
    end
end

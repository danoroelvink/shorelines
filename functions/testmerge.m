% x_mc=[500 500 500 700 900 900 900 700 500]
% y_mc=[500 700 900 900 900 700 500 500 500]

clear all;close all
alfa=360:-30:0;
x_mc=500+200*cosd(alfa);
y_mc=500+100*sind(alfa);
numb=[1:length(x_mc)];
if 1
    %clear all;close all
    figure(1);
    plot(x_mc,y_mc,'.-b','linewidth',2)
    text(x_mc,y_mc,num2str(numb'),'color','b','fontsize',16)
    axis([200 800 100 700]);
    hold on;
    [x,y]=select_polygon
    numb=[1:length(x)]
    text(x,y,num2str(numb'),'color','r','fontsize',16)
    
end
figure(2);
x_mc=[x_mc,nan,x];
y_mc=[y_mc,nan,y];
[ x_mc,y_mc ] = merge_coastlines_mc( x_mc,y_mc )
numb=[1:length(x_mc)];
plot(x_mc,y_mc,'k','linewidth',2)
text(x_mc,y_mc,num2str(numb'),'color','k','fontsize',16)
axis([200 800 100 700]);

% s(1)=0;
% for i=2:length(x);
%     s(i)=s(i-1)+hypot(x(i)-x(i-1),y(i)-y(i-1));
% end
% eps=.1;
%
% P=InterX([x;y]);
%
% plot(x,y,'.-b',P(1,:),P(2,:),'.k')
% iX=0;
% ind=[];
% for i=1:length(s)-1
%     for ip=1:size(P,2)
%         err=abs(s(i+1)-s(i)-hypot(P(1,ip)-x(i)  ,P(2,ip)-y(i)) ...
%                 -hypot(P(1,ip)-x(i+1),P(2,ip)-y(i+1)))
%         if err<eps
%             iX=iX+1;
%             ind(iX)=i;
%         end
%     end
% end
% if isempty(ind)
%     xnew=x;
%     ynew=y;
% else
% xnew=[x(1:ind(1)),x(ind(4)+1:end),nan,x(ind(2)+1:ind(3)),x(ind(2)+1),nan];
% ynew=[y(1:ind(1)),y(ind(4)+1:end),nan,y(ind(2)+1:ind(3)),y(ind(2)+1),nan];
% end
% plot(xnew,ynew,'k','linewidth',2)


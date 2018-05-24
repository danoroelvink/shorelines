function [ xnew,ynew ] = merge_coastlines( x,y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
eps=.1;
s(1)=0;
yesplot=0;
for i=2:length(x);
    if isnan(x(i))
        s(i)=nan;
    elseif isnan(x(i-1))
        s(i)=0;
    else
        s(i)=s(i-1)+hypot(x(i)-x(i-1),y(i)-y(i-1));
    end
end

P=InterX([x;y]);
if yesplot
    figure(3);
    plot(x,y,'.-b',P(1,:),P(2,:),'ok');
    hold on
    num=[1:length(x)];
    for i=1:length(x);
        text(x(i),y(i),num2str(num(i)));
    end
end
iX=0;
ind=[];
for i=1:length(s)-1
    if ~isnan(s(i))&&~isnan(s(i+1))
        for ip=1:size(P,2)
            err=abs(s(i+1)-s(i)-hypot(P(1,ip)-x(i)  ,P(2,ip)-y(i)) ...
                -hypot(P(1,ip)-x(i+1),P(2,ip)-y(i+1)));
            if err<eps
                iX=iX+1;
                ind(iX)=i;
            end
        end
    end
end
ind;
if isempty(ind) || length(ind)<2
    xnew=x;
    ynew=y;
elseif length(ind)<4
    xnew=[x(1:ind(1)),x(ind(2)+1:end),nan,x(ind(1)+1:ind(2)),x(ind(1)+1)];
    ynew=[y(1:ind(1)),y(ind(2)+1:end),nan,y(ind(1)+1:ind(2)),y(ind(1)+1)];
else
    xnew=[x(1:ind(1)),x(ind(4)+1:end),nan,x(ind(2)+1:ind(3)),x(ind(2)+1)];
    ynew=[y(1:ind(1)),y(ind(4)+1:end),nan,y(ind(2)+1:ind(3)),y(ind(2)+1)];
end
if yesplot
    plot(xnew,ynew,'k','linewidth',2)
    hold off
    if length(ind)>1
        ind
    end
end

end


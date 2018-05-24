function [ s,x,y ] = make_sgrid( x,y,ds0 )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
s=zeros(size(x));
s(1)=0;
for i=2:length(x)
    s(i)=s(i-1)+hypot(x(i)-x(i-1),y(i)-y(i-1));
end
snew=s;
%disp(['total length = ',num2str(s(end))]);
i=2;
while i<=length(snew)
    ds=snew(i)-snew(i-1);
    if ds<ds0/2
        %throw out point i
        snew=[snew(1:i-1),snew(i+1:end)];
    elseif ds>ds0*2
        %insert point i
        snew=[snew(1:i-1),.5*(snew(i-1)+snew(i)),snew(i:end)];
        i=i+1;
    else
        i=i+1;
    end
end
snew(2:end-1)=.25*snew(1:end-2)+.5*snew(2:end-1)+.25*snew(3:end);
x=interp1(s,x,snew);
y=interp1(s,y,snew);
s=snew;

end


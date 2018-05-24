function [s,x,y,x_mc,y_mc,ar]=make_sgrid_mc(x_mc,y_mc,ds0,i_mc,smoothfac)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
ar=0;
%% find number of sections
[ x,y,~,i1,i2 ] = get_one_polygon( x_mc,y_mc,i_mc );
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
snew(2:end-1)=smoothfac*snew(1:end-2)+(1.-2*smoothfac)*snew(2:end-1)+smoothfac*snew(3:end);
x=interp1(s,x,snew);
y=interp1(s,y,snew);
%ar=round(polyarea(x,y));
%disp(['section',num2str(i_mc),' area ',num2str(ar)]);
s=snew;
%% insert x,y back into x_mc,y_mc
x_mc=[x_mc(1:i1-1),x,x_mc(i2+1:end)];
y_mc=[y_mc(1:i1-1),y,y_mc(i2+1:end)];

end


function [x_mc,y_mc]=initialize_grid(x_mc,y_mc,ds0);
%% find number of sections
nans=find(isnan(x_mc));
n_mc=length(nans)+1;
yesplot=false;
if yesplot
   figure(11)
end
for i_mc=1:n_mc
    [ x,y,n_mc,i1,i2 ] = get_one_polygon( x_mc,y_mc,i_mc );   
    if yesplot
       plot(x,y); hold on
    end
    x0=x;
    y0=y;
    s0=zeros(size(x));
    for i=2:length(x)
        s0(i)=s0(i-1)+hypot(x(i)-x(i-1),y(i)-y(i-1));
    end
    ns=ceil(s0(end)/ds0);
    ds=s0(end)/ns;
    s=[0:ds:s0(end)];
    x=interp1(s0,x0,s);
    y=interp1(s0,y0,s);
    %% insert x,y back into x_mc,y_mc
    x_mc=[x_mc(1:i1-1),x,x_mc(i2+1:end)];
    y_mc=[y_mc(1:i1-1),y,y_mc(i2+1:end)];
    
end
if yesplot
   plot(x_mc,y_mc,'o')
end
end

function sedero_bargraph(x,y,dn,exagfac)

x=x/1000;
y=y/1000;
dn=dn/1000;
dX=zeros(size(x));
dY=zeros(size(x));
Hyp=zeros(size(x));
for i=2:length(x)-1;
    dX(i)=x(i+1)-x(i-1);
    dY(i)=y(i+1)-y(i-1);
    Hyp(i)=hypot(dX(i),dY(i));
end
plot(x,y,'k','linewidth',2)
axis equal
for i=2:length(x)-1
    xi(1)=.5*(x(i-1)+x(i));
    yi(1)=.5*(y(i-1)+y(i));
    xi(2)=.5*(x(i)+x(i+1));
    yi(2)=.5*(y(i)+y(i+1));
    cosa=dX(i)/Hyp(i);
    sina=dY(i)/Hyp(i);
    xi(3)=xi(2)-sina*dn(i)*exagfac;
    yi(3)=yi(2)+cosa*dn(i)*exagfac;
    xi(4)=xi(1)-sina*dn(i)*exagfac;
    yi(4)=yi(1)+cosa*dn(i)*exagfac;
    xi(5)=xi(1);
    yi(5)=yi(1);
    if dn(i)>0
        col='g'
    else
        col='r'
    end
    patch(xi,yi,col)
end
xlabel ('Easting (km)')
ylabel ('Northing (km)')

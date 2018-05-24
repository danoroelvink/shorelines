function [S]=update_shoreline(S);
distmax=abs(S.seamin/S.seaslope);
dntide_mc=S.x_mc*0;
sedero=S.zg-S.zg0;
yesplot=false;
if yesplot
    figure(12);clf;
end
xg=S.xg;
yg=S.yg;
zg=S.zg;
for i=1:size(xg,2);
    for j=1:size(xg,1);
        [dist(j,i),ip]=dist_to_polyline(S.x_mc,S.y_mc,xg(j,i),yg(j,i),10000);
        if dist(j,i)<distmax
            dV=(sedero(j,i))*S.dsdn(j,i);
            dntide_mc(ip)=dntide_mc(ip)+dV/S.ds0/S.d/S.nt;
            plot([S.x_mc(ip),xg(j,i)],[S.y_mc(ip),yg(j,i)],'m');
            hold on
        end
    end
end
n_mc=length(find(isnan(S.x_mc)))+1;
x_mc_old=S.x_mc;
y_mc_old=S.y_mc;
for i_mc=1:n_mc
    x=get_one_polyvar(S.x_mc,i_mc);
    y=get_one_polyvar(S.y_mc,i_mc);
    dntide=get_one_polyvar(dntide_mc,i_mc);
    
    n=length(x)-1
    %% Cyclic or not ?
    cyclic = hypot(x(end)-x(1),y(end)-y(1))<S.ds0;
    for i=1:n
        if cyclic
            im1=mod2(i-1,n);
            ip1=mod2(i+1,n);
        else
            im1=max(i-1,1);
            ip1=min(i+1,n+1);
        end
        dn=dntide(i);
        dx(i)=-dn*(y(ip1)-y(im1))/hypot(x(ip1)-x(im1),y(ip1)-y(im1));
        dy(i)= dn*(x(ip1)-x(im1))/hypot(x(ip1)-x(im1),y(ip1)-y(im1));
    end
    for i=1:n
        x(i)=x(i)+dx(i);
        y(i)=y(i)+dy(i);
    end
    if cyclic
        x(n+1)=x(1);
        y(n+1)=y(1);
    end
    [S.x_mc,S.y_mc]=insert_section(x,y,S.x_mc,S.y_mc,i_mc);
    
end

S.dntide=dntide_mc;
pcolor(xg,yg,sedero);
cmax=max(max(abs(sedero)));
caxis([-cmax cmax])
colorbar;
shading flat;
hold on
plot(x_mc_old,y_mc_old,S.x_mc,S.y_mc,'linewidth',2);
axis equal
if isfield(S,'ifig')
    S.ifig=S.ifig+1;
else
    S.ifig=1000;
end
% fname=[S.outputdir,filesep,'sedero',num2str(S.ifig),'.jpg'];
% print('-djpeg',fname)

%pause

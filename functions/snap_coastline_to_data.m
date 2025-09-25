xm=[100,200,400,300];
ym=[20,30,-10,-20];
xd=[0,350,1000];
yd=[0,0,0];
ds0=10;
dsf=.001;

%% original coastline as clicked xd,yd
ds=hypot(diff(xd),diff(yd));
sd=[0,cumsum(ds)];
%% Subdivide coastline with resolution ds0
ns=ceil(sd(end)/ds0);
ds=sd(end)/ns;
sc=[0:ds:sd(end)];
xc=interp1(sd,xd,sc);
yc=interp1(sd,yd,sc);
%% Subdivide coastline with high resolution dsf
ns=ceil(sd(end)/dsf);
ds=sd(end)/ns;
sf=[0:ds:sd(end)];
xf=interp1(sd,xd,sf);
yf=interp1(sd,yd,sf);
%% Compute distances from fine coastline xf,yf to measured xm,ym
dist=hypot(xm-xf',ym-yf');
[distmin,ind]=min(dist);
dirc=atan2d(yf(ind+1)-yf(ind-1),xf(ind+1)-xf(ind-1));
dirm=atan2d(ym-yf(ind-1),xm-xf(ind-1));
dmin=distmin.*sign(sind(dirm-dirc));
% Sort by ascending s 
[ind,ind2]=sort(ind);
dmin=dmin(ind2);
corr=interp1(sf(ind),dmin,sc);
corr(isnan(corr))=0;
figure;plot(xc,yc,xm,ym,'o',xc,corr)
dirc=zeros(size(xc));
dirc(1)=atan2d(yc(2)-yc(1),xc(2)-xc(1));
dirc(end)=atan2d(yc(end)-yc(end-1),xc(end)-xc(end-1));
dirc(2:end-1)=atan2d(yc(3:end)-yc(1:end-2),xc(3:end)-xc(1:end-2));
dx= corr.*sind(dirc);
dy= corr.*cosd(dirc);
xccorr=xc+dx;
yccorr=yc+dy;
figure;plot(xc,yc,xm,ym,'o',xccorr,yccorr,xccorr,yccorr,'.')
axis equal


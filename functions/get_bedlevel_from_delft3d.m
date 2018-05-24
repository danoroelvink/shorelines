function [x,y,zb0,zb,dsdn]=get_bedlevel_from_delft3d(fname)
vs_use(fname);
times=vs_let('map-info-series','ITMAPC');
last=length(times);
zb0=-squeeze(vs_let('map-sed-series',{1},'DPS','quiet'));
zb=-squeeze(vs_let('map-sed-series',{last},'DPS','quiet'));
dsdn=-squeeze(vs_let('map-const','GSQS','quiet'));
sedero=zb-zb0;
x=squeeze(vs_let('map-const','XZ','quiet'));
y=squeeze(vs_let('map-const','YZ','quiet'));
y(x==0)=nan;
zb(x==0)=nan;
x(x==0)=nan;
if 1
figure;
pcolor(x,y,zb0)
shading interp;
axis equal;
colorbar
end
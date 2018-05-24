function [S]=update_bathy(S);

xg=S.xg;
yg=S.yg;
zg0=S.zg;
for i=1:size(xg,2);
    for j=1:size(xg,1);
        [dist(j,i),ip]=dist_to_polyline(S.x_mc,S.y_mc,xg(j,i),yg(j,i),10000);
    end
end
IN=inpolygon(xg,yg,S.x_mc,S.y_mc);
zg=-dist*S.seaslope;
zg(zg<S.seamin)=zg0(zg<S.seamin);
zg(IN)=min(dist(IN)*S.landslope,S.landmax);
zg(isnan(zg))=999;
dps=-zg;
save(S.bathyname,'dps','-ascii')
save('x.dep','xg','-ascii')
S.zg=zg;

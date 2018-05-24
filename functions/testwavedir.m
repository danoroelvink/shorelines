addpath(genpath('..\..\..\functions\'))
scale=1000;
fname='xboutput.nc';
x=nc_varget(fname,'globalx');
y=nc_varget(fname,'globaly');
H=nc_varget(fname,'H',[1 0 0],[1 Inf Inf]);
thetamean=nc_varget(fname,'thetamean',[1 0 0],[1 Inf Inf]);
figure
pcolor(x,y,H);
shading interp;
axis equal;
colorbar;
hold on;
[x_mc,y_mc]=select_multi_polygon('k');
for i=1:length(x_mc);
    if ~isnan(x_mc(i));
        [row,col]=find_in_grid(x,y,x_mc(i),y_mc(i));
        mag=H(row,col);
        dir=thetamean(row,col);
        vecx=[x_mc(i) x_mc(i)-scale*mag*sin(dir*pi/180)];
        vecy=[y_mc(i) y_mc(i)-scale*mag*cos(dir*pi/180)];
        plot(vecx,vecy,'w',vecx(1),vecy(1),'ow','linewidth',2)
    end
end

fnames1={'stlouis_2004_back_bay_north_utm28n.ldb',...
    'stlouis_2004_back_bay_south_utm28n.ldb',...
    'stlouis_2004_north_spit_utm28n.ldb',...
    'stlouis_2004_south_spit_utm28n.ldb'};
fnames2={'stlouis_2016_back_bay_north_utm28n.ldb',...
    'stlouis_2016_back_bay_south_utm28n.ldb',...
    'stlouis_2016_north_spit_utm28n.ldb',...
    'stlouis_2016_south_spit_utm28n.ldb'};
figure(1)
subplot(121)
x_mc=[];
y_mc=[];
for i=1:length(fnames1)
    xy=load(fnames1{i});
    if i~=3
        x=xy(:,1);
        y=xy(:,2);
    else
        x=xy(end:-1:1,1);
        y=xy(end:-1:1,2);
    end
    x_mc=[x_mc;nan;x];
    y_mc=[y_mc;nan;y];
    plot(x,y,x(1),y(1),'o','linewidth',2)
    axis equal
    hold on;
end
out(:,1)=x_mc(2:end);
out(:,2)=y_mc(2:end);
save('2004.ldb','out','-ascii')
subplot(122)
x_mc=[];
y_mc=[];
for i=1:length(fnames2)
    xy=load(fnames2{i});
    if i~=3
        x=xy(:,1);
        y=xy(:,2);
    else
        x=xy(end:-1:1,1);
        y=xy(end:-1:1,2);
    end
    x_mc=[x_mc;nan;x];
    y_mc=[y_mc;nan;y];
    plot(x,y,x(1),y(1),'o','linewidth',2)
    axis equal
    hold on;
end
out=[];
out(:,1)=x_mc(2:end);
out(:,2)=y_mc(2:end);
save('2016.ldb','out','-ascii')


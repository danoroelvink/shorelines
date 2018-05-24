%% Analyze results of shorelines simulation
%addpath(genpath('d:\data\dano\ShorelineS\'))
data_end='..\ijmuiden\xy2007.txt';
modelresults='..\Results.mat';
obstitle='Observed coastline changes 1967-2007';
modtitle='Modelled coastline changes 1967-2007';

ylims=[480 510];
width=1500;
exagfac=3;

load(modelresults);
nt=length(S.xyt);
x_mc0=S.xyt(1).x_mc;
y_mc0=S.xyt(1).y_mc;
%% Construct normal
figure;
plot(x_mc0,y_mc0,'k.-')
axis equal
hold on
x_n=zeros(length(x_mc0),2)
y_n=zeros(length(x_mc0),2)
dX=zeros(size(x_mc0));
dY=zeros(size(x_mc0));
Hyp=zeros(size(x_mc0));
for i=2:length(x_mc0)-1;
    dX(i)=x_mc0(i+1)-x_mc0(i-1);
    dY(i)=y_mc0(i+1)-y_mc0(i-1);
    Hyp(i)=hypot(dX(i),dY(i));
    dx=-width*dY(i)/Hyp(i);
    dy= width*dX(i)/Hyp(i);
    x_n(i,:) = [x_mc0(i)-dx,x_mc0(i)+dx];
    y_n(i,:) = [y_mc0(i)-dy,y_mc0(i)+dy];
end
ds0=Hyp(2)/2;
x_n(x_n==0)=nan;
y_n(y_n==0)=nan;
plot(x_n',y_n','k',x_n(:,1)',y_n(:,1)','o');

x_mc=S.xyt(end).x_mc;
y_mc=S.xyt(end).y_mc;
for i=2:length(x_mc0)-1
    L1=[x_mc;y_mc];
    L2=[squeeze(x_n(i,:));squeeze(y_n(i,:))];
    P{i}=InterX(L1,L2);
end
xend=x_mc0;
yend=y_mc0;

for i=2:length(x_mc0)-1
    p=P{i};
    if ~isempty(p)
        xend(i)=p(1);
        yend(i)=p(2);
    else
        xend(i)=nan;
        yend(i)=nan;
    end
end
plot(xend,yend,'r')
dn=hypot(xend-x_n(:,1)',yend-y_n(:,1)')-width;

xy=load(data_end);
x_mc=xy(:,1)';
y_mc=xy(:,2)';
for i=2:length(x_mc0)-1
    L1=[x_mc;y_mc];
    L2=[squeeze(x_n(i,:));squeeze(y_n(i,:))];
    P{i}=InterX(L1,L2);
end
xendm=x_mc0;
yendm=y_mc0;

for i=2:length(x_mc0)-1
    p=P{i};
    if ~isempty(p)
        xendm(i)=p(1);
        yendm(i)=p(2);
    else
        xendm(i)=nan;
        yendm(i)=nan;
    end
end
dnm=hypot(xendm-x_n(:,1)',yendm-y_n(:,1)')-width;
figure;
plot(dn,'r')
hold on;
plot(dnm,'b')
rmserr=std(dn(~isnan(dnm))-dnm(~isnan(dnm)));
bias=mean(dn(~isnan(dnm))-dnm(~isnan(dnm)));
text(20,500,['rms error = ',num2str(rmserr),'  bias = ',num2str(bias)]);
figure
subplot(121)
sedero_bargraph(x_mc0,y_mc0,dnm,exagfac)
ylim(ylims)
title(obstitle)
subplot(122)
sedero_bargraph(x_mc0,y_mc0,dn,exagfac)
ylim(ylims)
title(modtitle)

print('-dpng','-r300','bargraphs.png')



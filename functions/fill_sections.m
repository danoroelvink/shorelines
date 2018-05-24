function [xmax,ymax,xmin,ymin,nmax ] = fill_sections(x_mc,y_mc,ld,it,xmax,ymax,xmin,ymin,nmax)
%function [nmax] = fill_sections(x_mc,y_mc,ld,it,nmax)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nans=find(isnan(x_mc));
n_mc=length(nans)+1;

for i_mc=1:n_mc
    nans=find(isnan(x_mc));
    if isempty(nans)
        i1=1;
        i2=length(x_mc);
    else
        if i_mc==1
            i1=1;
            i2=nans(i_mc)-1;
        elseif i_mc==n_mc;
            i1=nans(i_mc-1)+1;
            i2=length(x_mc);
        else
            i1=nans(i_mc-1)+1;
            i2=nans(i_mc)-1;
        end
    end
    x=x_mc(i1:i2);
    y=y_mc(i1:i2);
    phi=[];
    for i=1:length(x)-1
        phi(i)=atan2((y(i+1)-y(i)),(x(i+1)-x(i)));
    end

    % MAKE SURE TO FILL LAND AND LAKES DIFFERENTLY
    dphi=[];
    for i=1:length(x)-2
        dphi(i)=phi(i+1)-phi(i);
        if dphi(i)>pi
            dphi(i)=dphi(i)-2*pi;
        end
        if dphi(i)<-pi
            dphi(i)=dphi(i)+2*pi;
        end
    end
    %right=sum(dphi)<(1.8*pi);                                              % theoretically a lake should have : sum(dphi)==2*pi (but in practice not exactly)
    right=sum(dphi)<0;                                              % theoretically a lake should have : sum(dphi)==2*pi (but in practice not exactly)

    notclosed=1;   % <- notclosed option recently added for Ameland (check?)
    if (x(1)-x(end)).^2+(y(1)-y(end)).^2 == 0
    notclosed=0;
    end
    
    if 0
    % CONSTRUCT A LAND POLYGON WHICH DOES NOT CROSS THE COASTLINE
    if (right || i_mc==1) && notclosed  % <- notclosed option recently added for Ameland (check?)
        dx=diff(x);dx=[dx(1),(dx(1:end-1)+dx(2:end))/2,dx(end)];
        dy=diff(y);dy=[dy(1),(dy(1:end-1)+dy(2:end))/2,dy(end)];
        l2=(dx.^2+dy.^2).^0.5;
        xclose = [x(end)+dy(end)/l2(end),x(1)+dy(1)/l2(1)];
        yclose = [y(end)-dx(end)/l2(end),y(1)-dx(1)/l2(1)];
      
        
        %plot(x_mc,y_mc,dx,dy,xclose,yclose,'*')
        % look for points left and right of line
        RC = diff(yclose)/diff(xclose);
        if diff(xclose)==0
           RC = diff(yclose)/1e-5;
        end
        B = yclose(1)-RC*xclose(1);
        ID = find(y<RC*x+B);                                       
        %ID = find(y>RC*x+B);
        %[xcr,ycr,x1n,y1n,x2n,y2n,ID]=findCrossings2(x,y,xclose,yclose);       % <- more sophisticated approach
        l3=100;
        xclose = [xclose(1),fliplr(x(ID)+l3*dy(ID)./l2(ID)),xclose(2)];
        yclose = [yclose(1),fliplr(y(ID)-l3*dx(ID)./l2(ID)),yclose(2)];
        x = [x,xclose];
        y = [y,yclose];
    end
    end
%% For straight line add land fill behind shoreline
if notclosed
[x,y,xmax,ymax,xmin,ymin,nmax] = insert_land_fill(x,y,ld,it,xmax,ymax,xmin,ymin,i_mc,n_mc,nmax);
%[x,y,nmax] = insert_land_fill(x,y,ld,it,i_mc,n_mc,nmax);
end
%%
    % fill points
    if right || i_mc==1
        fill(x,y,[.8 .8 0]);
    else
        fill(x,y,[.9 .9 1]);
    end

    hold on
end
hold off


%end


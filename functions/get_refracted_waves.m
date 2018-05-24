function [H phiw_cd]=get_refracted_waves(x,y,surf_width,xg,yg,Hg,phiwg)
scale=1000;
yesplot=1;
H=[];
phiw_cd=[];
if yesplot
    figure(15)
    pcolor(xg,yg,Hg);
    shading flat;
    axis equal;
    colorbar;
    hold on;
end
%% create line at depth of closure x_cd, y_cd
n=length(x)-1;
for i=1:n
    dX=x(i+1)-x(i);
    dY=y(i+1)-y(i);
    Hyp=hypot(dX,dY);
    dx=-surf_width*dY/Hyp;
    dy= surf_width*dX/Hyp;
    x_cd(i)=0.5*(x(i)+x(i+1))+dx;
    y_cd(i)=0.5*(y(i)+y(i+1))+dy;
end
for i=1:n;
    if ~isnan(x);  
        [row,col]=find_in_grid(xg,yg,x_cd(i),y_cd(i));
        if (size(row,1)>1)
            continue%/break to check
        end
%         try
%             disp(['i ',num2str(i),' row ',num2str(row),' col ',num2str(col)])
%         catch
%             figure(95)
%             plot(x(1:154),y(1:154),'-o',x_cd(1:154),y_cd(1:154),'-o',xg,yg,'k',xg',yg','k')
%             x(154)
%             y(154)
%             x_cd(154)
%             y_cd(154)
%             i
%             row(1:10)
%             col(1:10)
%         end
        H(i)=Hg(row,col);
        phiw_cd(i)=phiwg(row,col);%*pi/180;
        if yesplot
            vecx=[x_cd(i) x_cd(i)-scale*H(i)*sin(phiw_cd(i))];
            vecy=[y_cd(i) y_cd(i)-scale*H(i)*cos(phiw_cd(i))];
            plot(vecx,vecy,'k',vecx(1),vecy(1),'ok','linewidth',2)
        end
    end
end


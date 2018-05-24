function [ x_mc,y_mc ] = merge_coastlines_mc(x_mc,y_mc)
yesplot=0;
eps=1;
[ x,y,n_mc,i1,i2 ] = get_one_polygon( x_mc,y_mc,1 );
for i_mc=1:n_mc
    for j_mc=1:n_mc
        if i_mc~=j_mc
            [ xi,yi,n_mc,i1i,i2i ] = get_one_polygon( x_mc,y_mc,i_mc );
            [ xj,yj,n_mc,i1j,i2j ] = get_one_polygon( x_mc,y_mc,j_mc );
            if length(xi)>2&&length(xj)>2
                si=cumdist(xi,yi);
                sj=cumdist(xj,yj);
                P=InterX([xi;yi],[xj;yj]);
                iX=0;
                indi=[];
                for i=1:length(si)-1
                    for ip=1:size(P,2)
                        err=abs(si(i+1)-si(i)-hypot(P(1,ip)-xi(i)  ,P(2,ip)-yi(i)) ...
                            -hypot(P(1,ip)-xi(i+1),P(2,ip)-yi(i+1)));
                        if err<eps
                            iX=iX+1;
                            indi(iX)=i;
                        end
                    end
                end
                if yesplot
                    figure(3)
                    plot(xi,yi,'.-b',xj,yj,'.-r',P(1,:),P(2,:),'ok');
                    hold on
                    num=[1:length(xi)];
                    for i=1:length(xi);
                        text(xi(i),yi(i),num2str(num(i)));
                    end
                    num=[1:length(xj)];
                    for i=1:length(xj);
                        text(xj(i),yj(i),num2str(num(i)));
                    end
                    hold off
                end
                
                iX=0;
                indj=[];
                for i=1:length(sj)-1
                    for ip=1:size(P,2)
                        err=abs(sj(i+1)-sj(i)-hypot(P(1,ip)-xj(i)  ,P(2,ip)-yj(i)) ...
                            -hypot(P(1,ip)-xj(i+1),P(2,ip)-yj(i+1)));
                        if err<eps
                            iX=iX+1;
                            indj(iX)=i;
                        end
                    end
                end
                if isempty(indi)
                    xnewi=xi;
                    ynewi=yi;
                    xnewj=xj;
                    ynewj=yj;
                elseif length(indi)>=2
                    if indi(1)>=1&&indj(1)>=1
                        xnewi=[xi(1:indi(1)),P(1,1),xj(indj(2)+1:end),xj(1:indj(1)),P(1,2),xi(indi(2)+1:end)];
                        ynewi=[yi(1:indi(1)),P(2,1),yj(indj(2)+1:end),yj(1:indj(1)),P(2,2),yi(indi(2)+1:end)];
                        xnewj=[xi(1:indi(1)),P(1,2),xj(indj(2)+1:end),xj(1:indj(1)),P(1,1),xi(indi(2)+1:end)];
                        ynewj=[yi(1:indi(1)),P(2,2),yj(indj(2)+1:end),yj(1:indj(1)),P(2,1),yi(indi(2)+1:end)];
                        snewi=cumdist(xnewi,ynewi);
                        snewj=cumdist(xnewj,ynewj);
                        if snewj(end)<snewi(end)
                            xnewi=xnewj;
                            ynewi=ynewj;
                        end
                    elseif indi(1)==1
                        xnewi=[P(1,1),xi(2:indi(2)),P(1,2),xj(indj(2)+1:end),xj(1:indj(1)),P(1,1)];
                        ynewi=[P(2,1),yi(2:indi(2)),P(2,2),yj(indj(2)+1:end),yj(1:indj(1)),P(2,1)];
                        xnewj=[P(1,2),xi(2:indi(2)),P(1,1),xj(indj(2)+1:end),xj(1:indj(1)),P(1,2)];
                        ynewj=[P(2,2),yi(2:indi(2)),P(2,1),yj(indj(2)+1:end),yj(1:indj(1)),P(2,2)];
                        snewi=cumdist(xnewi,ynewi);
                        snewj=cumdist(xnewj,ynewj);
                        if snewj(end)<snewi(end)
                            xnewi=xnewj;
                            ynewi=ynewj;
                        end
                    elseif indj(1)==1
                        xnewi=[P(1,1),xj(2:indj(2)),P(1,2),xi(indi(2)+1:end),xi(1:indi(1)),P(1,1)];
                        ynewi=[P(2,1),yj(2:indj(2)),P(2,2),yi(indi(2)+1:end),yi(1:indi(1)),P(2,1)];
                        xnewj=[P(1,2),xj(2:indj(2)),P(1,1),xi(indi(2)+1:end),xi(1:indi(1)),P(1,2)];
                        ynewj=[P(2,2),yj(2:indj(2)),P(2,1),yi(indi(2)+1:end),yi(1:indi(1)),P(2,2)];
                        snewi=cumdist(xnewi,ynewi);
                        snewj=cumdist(xnewj,ynewj);
                        if snewj(end)<snewi(end)
                            xnewi=xnewj;
                            ynewi=ynewj;
                        end
                        
                    end
                    if yesplot
                        hold on
                        plot(xnewi,ynewi,'k','linewidth',2)
                        hold off
                        
                        indi;
                        indj;
                        pause
                        
                    end
                    xnewj=[];
                    ynewj=[];
                else
                    xnewi=[];
                    ynewi=[];
                    xnewj=[];
                    ynewj=[];
                end
                [x_mc,y_mc]=insert_section(xnewi,ynewi,x_mc,y_mc,i_mc);
                [x_mc,y_mc]=insert_section(xnewj,ynewj,x_mc,y_mc,j_mc);
            end
        end
    end
    
end


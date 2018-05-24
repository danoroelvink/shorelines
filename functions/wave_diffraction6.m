function [QS,xp,yp]=wave_diffraction2(QS,x,y,S,hsbr,sphibr,phiw,x_mc,y_mc,x_hard,y_hard,n,hbr)
%% wave diffraction following Roelvink approach for angles, Kamphuis for Kd
if ~isempty(x_hard)
    
    A=0.1; % for Dean profile
    n_hard=length(find(isnan(x_hard)))+1;
    Kdrs(1:n_hard,1:n)=1;
    % phibr_cor(1:n_hard,1:n)=1;
    % krt(1:n_hard,1:n)=1;
    kdt(1:n_hard,1:n)=1;
    % ks(1:n_hard,1:n)=1;
    kd(1,1:2)=1;
    
    
    xS=.5*(x(1:n)+x(2:n+1));
    yS=.5*(y(1:n)+y(2:n+1));
    len=5*hypot(max(x_mc)-min(x_mc),max(y_mc)-min(y_mc));
    for i=1:n
        
        
        xw=[xS(i)-len*sin(phiw),xS(i)+len*sin(phiw)]; % should it be x or x_mc or xS
        yw=[yS(i)-len*cos(phiw),yS(i)+len*cos(phiw)];
        
        
        n_hard=length(find(isnan(x_hard)))+1;
        
        %% check if the point in the structure shadow
        for i_hard=1:n_hard
            Ps1=zeros(2,2);

            [ x_struc,y_struc,n_hard,~,~ ] = get_one_polygon( x_hard,y_hard,i_hard );

            %catch structure tip points
            xtip(1)=x_struc(1);
            ytip(1)=y_struc(1);
            xtip(2)=x_struc(end);
            ytip(2)=y_struc(end);
             
            %tip depth & wave length
            dtip =A*(mean(ytip)).^(2/3); % should be adjusted for other orientations
            
            %check if the tip in the breaking / broken zone
            if dtip > hbr(i)
                %calculate the wave height at the tip
                [~,c0]=GUO2002(S.tper,S.ddeep);
                [ktip,ctip]=GUO2002(S.tper,dtip);
                phi_tip=asin((ctip/c0).*sin(phiw));
                cg0=c0/2;
                cgtip=ctip*(0.5+ktip*dtip/sinh(2*ktip*dtip));
                hstip=S.Hso*sqrt(cg0*cos(phiw)/cgtip/cos(phi_tip));
                
            else
                hstip=hbr(i)*S.gamma;
            end
            
            Lo =S.g*S.tper.^2./2./pi;
            Ltip = Lo.*sqrt(tanh(4*pi.^2.*dtip./S.tper.^2./S.g));
            
            G=2.5*Ltip; %width of the transition zone at shoreline
            
            
            disd=(hbr./A).^(3/2);   %distance between the shoreline and the breaking depth

            
            %% the point (P) where we calculate the coefficients (xp,yp)
            % approach #3
            dX=x(i+1)-x(i);
            dY=y(i+1)-y(i);
            Hyp=hypot(dX,dY);
            dx=-disd*dY/Hyp;
            dy= disd*dX/Hyp;
            xp(i)=0.5*(x(i)+x(i+1))+dx(i);
            yp(i)=0.5*(y(i)+y(i+1))+dy(i);
            
            if yp(i)==min(ytip);
                yp(i)=yp(i)-0.5;
            end
            
            % interpolate 2 tip projection points
            for itip=1:2
                x_meged(1)=xtip(itip);
                y_meged(1)=ytip(itip);
                x_meged(2)=xtip(itip)-len*sin(phiw);
                y_meged(2)=ytip(itip)-len*cos(phiw);
                temp=InterX([xS;yS],[x_meged;y_meged]);
                
                % to tackel the problem of multi intersect points in case of salient / tombolo
                if size(temp,2)>1
                    for ifp=1:size(temp,2)
                        h(ifp)=hypot((xtip(1)-temp(1,ifp)),(ytip(1)-temp(2,ifp)));
                    end
                    k=find(min(h));
                else
                    k=1;
                end
                
                Ps1(:,itip)=temp(:,k); %Ps1(:,itip)=unique(ceil(temp)','rows')';
            end
            
            % extend the shadow area
            Ps2=Ps1;
            Ps2(1,1)= Ps1(1,1)-G; %should be adjusted to fit all oreintations
            Ps2(1,2)= Ps1(1,2)+G;
            
            
            % check if the wave ray intersect with the projected(shadowed) line
            P=InterX([xw;yw],[Ps2(1,:);Ps2(2,:)]);
            
            
            % Here we go!!!
            shadowed(i)=size(P,2)>0;
            if size(P,2)>0
                
                %% Diffraction Coefficient
                % calculate shadow line angle & theta
                for itip=1:2
                    alpham2=atan2(real((Ps2(2,itip)-ytip(itip))),real(( Ps2(1,itip)-xtip(itip))))*180/pi;
                    alpham=atan2(real((Ps1(2,itip)-ytip(itip))),real(( Ps1(1,itip)-xtip(itip))))*180/pi;
                    alphas=atan2(yp(i)-(ytip(itip)),(xp(i)-xtip(itip)))*180/pi;
                    delta=(alpham-alpham2);
                    
                    if itip == 1
                        theta=-(alphas-alpham);
                        alpahs_cor=alphas+delta;
                        if  alpahs_cor > 0
                            alpahs_cor =0
                        end
                    elseif itip == 2
                        theta=(alphas-alpham);
                        alpahs_cor=alphas+delta;
                        if  alpahs_cor < -180
                            alpahs_cor =-180
                        end
                    end
                    
                    if theta >= 0 && theta < abs(delta)
                        if itip ==1
                            alpahs_cor=alpham-theta+delta;
                        else itip==2
                            alpahs_cor=alpham+delta+theta;
                        end
                    end
                    
                    if theta > abs(delta)
                        alpahs_cor=-phiw*180/pi+270;
                    end
                    
                    
                    alpahs_x(i,itip)=pi+alpahs_cor*pi/180; %wave angle to x-axis
                    alpahs_N(i,itip)=1.5*pi-(alpahs_cor*pi/180); %wave angle to North
                    
                    % refraction coefficient
                    kr(itip)=real(sqrt(cos(alpahs_N(i,itip))/cos(phiw)));
                    
                    % calculate diffraction coeffcient kd from each tip
                    if  theta <= 0 && theta >= -90
                        kd(itip)=abs(0.69+0.008*theta*(1));
                    elseif theta > 0 && theta <= 40
                        kd(itip)=0.71+0.37*sin(theta*(pi/180));
                    elseif theta > 40 && theta <=90
                        kd(itip)=0.83+0.17*sin(theta*(pi/180));
                    elseif theta < -90
                        kd(itip)=0;
                    else
                        kd(itip)=0;
                    end
                end
                
                
                % combining kd from both tips
                kdt(i_hard,i)=sqrt(kd(1)^2+kd(2)^2+2*kd(1)*kd(2)*0);
                
                
                %% combining wave angles from both tips in the shadow zone
                hzl_cor=((kd(1)*hstip)^2)*cos(alpahs_x(i,1))+((kd(2)*hstip)^2)*cos(alpahs_x(i,2));
                vl_cor=((kd(1)*hstip)^2)*sin(alpahs_x(i,1))+((kd(2)*hstip)^2)*sin(alpahs_x(i,2));
                
                
                phicor_N=-(atan2(vl_cor,hzl_cor)+0.5*pi)+pi; %from North
                phic(i)=2*pi-atan2(y(i+1)-y(i),x(i+1)-x(i)); %coastline oreintation
                phicor(i)=atan2(sin(phic(i)-phicor_N),cos(phic(i)-phicor_N)); %normal to shoreline
                sphibr(i)=phicor(i);
                
                
                
                tombolo(i)=(mean(ytip)-y(i)) < 1;

                
                %% total Coeffcients
                Kdrs(i_hard,i)=kdt(i_hard,i)*1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% new conepct %%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%
        end
    end

   
    %% Combine multi structure effect
    if   n_hard > 1
        for i_hard=1:n_hard-1
            Kdrs(n_hard,:)= Kdrs(n_hard,:).*Kdrs(n_hard-i_hard,:);
        end
    end
    
    %% Calculate the long sore sediment transport (Overwrite the upwind !!change)
    for i=1:n
        im1=max(i-1,1);
        ip1=min(i+1,n);
        

    if shadowed(i)==0
        hs(i)=(hstip*Kdrs(n_hard,i));
    else
        hs(i)=hbr(i)*S.gamma;
    end

end

for i=1:n
    im1=max(i-1,1);
    ip1=min(i+1,n);
    
    QSkampmass(i)=2.33 * S.rhos/(S.rhos-S.rhow) .* S.tper.^1.5 .* S.tanbeta.^0.75 .* S.d50.^-0.25 .* (hs(i)).^2 .* ( (abs(sin(2*sphibr(i))).^0.6.*sign(sphibr(i)))-(2/S.tanbeta)*cos(sphibr(i)).*(hs(ip1)-hs(im1))/( hypot(yp(ip1)-yp(im1),(xp(ip1)-xp(im1)))));  %4*S.ds0
    
    QS(i) = 365*24*60*60*(QSkampmass(i) /(S.rhos-S.rhow)) /(1.0-S.porosity);
    
    QS(tombolo)=0;
end
end


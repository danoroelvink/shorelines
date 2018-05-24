function  [QS,xp,yp]=wave_diffraction4(QS,x,y,S,hsbr,sphibr,phiw,x_mc,y_mc,x_hard,y_hard,n,hbr);
%% wave diffraction following (hurst) approach
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
for i=1:n  %think again
    
    
    xw=[xS(i)-len*sin(phiw),xS(i)+len*sin(phiw)]; % should it be x or x_mc or xS
    yw=[yS(i)-len*cos(phiw),yS(i)+len*cos(phiw)];
    
    
    n_hard=length(find(isnan(x_hard)))+1;
    
    %% check if the point in the structure shadow
    for i_hard=1:n_hard
        Ps1=zeros(2,2);
        %         x_meged=zeros(1,2);
        %         y_meged=zeros(1,2);
        
        [ x_struc,y_struc,n_hard,~,~ ] = get_one_polygon( x_hard,y_hard,i_hard );
        
        
        
        %catch structure tip points
        xtip(1)=x_struc(1);
        ytip(1)=y_struc(1);
        xtip(2)=x_struc(end);
        ytip(2)=y_struc(end);
        
        
        %tip depth & wave length
        dtip =A*(mean(ytip)-0).^(2/3); % should be adjusted for other orientations
       %dtip=(mean(ytip)-y(i)).*S.tanbeta;
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
%             disd=(hbr./S.tanbeta);
           
 %% the point (P) where we calculate the coefficients (xp,yp)
               
                % approach #1

%                         xp(i)=x(i); %should be adjusted to fit all orientations
%                         yp(i)=disd(i)+y(i);
                % approach #2 
                
%                  xp(i)=x(i); 
%                  yp(i)=0;
%                  yp(i)=(mean(ytip)-1.5*Ltip);
%                  yp(i)=(mean(ytip)-S.ypd);
               
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
                     theta2=-(alphas-alpham2);
                     alpahs_cor=alpham-1.5*theta2;
                elseif itip == 2
                    theta=(alphas-alpham);
                    theta2=(alphas-alpham2);
                    alpahs_cor=alpham+1.5*theta2;
                end
%                 if theta >= 0 && theta < abs(delta)
%                     if itip ==1
%                         alpahs_cor=alpham-theta+delta;
%                     else itip==2
%                         alpahs_cor=alpham+delta+theta;
%                     end
%                 end
                
                if theta > abs(delta)
                    alpahs_cor=-phiw*180/pi+270;
                end
                
                %for hurst method 
                alpahs_x(i,itip)=pi+alpahs_cor*pi/180;  %wave angle to x-axis
                alpahs_N(i,itip)=1.5*pi-(alpahs_cor*pi/180); %wave angle to North
                % refraction coefficient
                kr(itip)=real(sqrt(cos(alpahs_N(i,itip))/cos(phiw)));
                
                % calculate diffraction coeffcient kd from each tip
                if  theta <= 0 && theta >= -90
%                     kd(itip)=0.71-0.0093*theta*(pi/180)+0.000025*(theta*(pi/180))^2;
                    kd(itip)=abs(0.69+0.008*theta*(1));
                    
%                     kd(itip)=0.5*cos(theta*(pi/180)); %hurst
                elseif theta > 0 && theta <= 40
                    kd(itip)=0.71+0.37*sin(theta*(pi/180));
                    %hurst
%                     ls=hypot(Ps2(2,1)-Ps2(2,2),Ps2(1,1)-Ps2(1,2));
%                     xs=hypot(Ps2(2,itip)-y(i),Ps2(1,itip)-x(i));
%                     kd(itip)=0.5*(1-sin(pi*xs/ls));
                elseif theta > 40 && theta <=90
                    kd(itip)=0.83+0.17*sin(theta*(pi/180));
%                     ls=hypot(Ps2(2,1)-Ps2(2,2),Ps2(1,1)-Ps2(1,2));
%                     xs=hypot(Ps2(2,itip)-y(i),Ps2(1,itip)-x(i));
%                     kd(itip)=0.5*(1-sin(pi*xs/ls));
                elseif theta < -90
                    kd(itip)=0;
                else
                    kd(itip)=0;
                end

            end
            
                     %             dp(i)=hbr(i);
            dp(i)=A*(hypot(yp(i)-y(i),xp(i)-x(i))).^(2/3);
            Lp(i) = Lo.*sqrt(tanh(4*pi.^2.*dp(i)./S.tper.^2./S.g));
            
            
            % Shoaling coefficient
            kp(i)=2*pi/Lp(i);
            np(i)=0.5*(1+(2*kp(i)*dp(i)/(sinh(2*kp(i)*dp(i)))));
            
            kI=2*pi/Ltip;
            nI=0.5*(1+(2*kI*dtip/(sinh(2*kI*dtip))));
            
            ks(i_hard,i)=sqrt((nI*Ltip)/(np(i)*Lp(i)));
            
            % combining kd from both tips
                kdt(i_hard,i)=sqrt(kd(1)^2+kd(2)^2+2*kd(1)*kd(2)*0);
            % combining kr from both tips 
%             krt(i_hard,i)=sqrt(kr(1)^2+kr(2)^2+2*kr(1)*kr(2)*0);

               %% combining wave angles from both tips in the shadow zone
            %for hurst method
            hzl_cor=((kd(1)*hstip)^2)*cos(alpahs_x(i,1))+((kd(2)*hstip)^2)*cos(alpahs_x(i,2));
            vl_cor=((kd(1)*hstip)^2)*sin(alpahs_x(i,1))+((kd(2)*hstip)^2)*sin(alpahs_x(i,2));
            
            
            phicor_N=-(atan2(vl_cor,hzl_cor)+0.5*pi)+pi; %from North
            phic(i)=2*pi-atan2(y(i+1)-y(i),x(i+1)-x(i)); %coastline oreintation
            phicor(i)=atan2(sin(phic(i)-phicor_N),cos(phic(i)-phicor_N)); %normal to shoreline

            sphibr(i)=phicor(i);

            tombolo(i)=(mean(ytip)-y(i)) < 1;
            krt(i_hard,i)=real(sqrt(cos(phicor_N)/cos(phiw)));
            %% total Coefficients
%                         Kdrs(i_hard,i)=krt(i_hard,i)*kdt(i_hard,i)*ks(i_hard,i);
            Kdrs(i_hard,i)=kdt(i_hard,i)*ks(i_hard,i);
       
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hbr2(i)=hstip*Kdrs(i_hard,i)/0.72;
disd2(i)=(hbr2(i)/A).^(3/2); 


                                    dX=x(i+1)-x(i);
                                    dY=y(i+1)-y(i);
                                    Hyp=hypot(dX,dY);
                                    dx(i)=-disd2(i)*dY/Hyp;
                                    dy(i)= disd2(i)*dX/Hyp;
                                    xp(i)=0.5*(x(i)+x(i+1))+dx(i);
                                    yp(i)=0.5*(y(i)+y(i+1))+dy(i);

%         xp(i)=x(i); %should be adjusted to fit all orientations
if yp(i)==min(ytip);
yp(i)=yp(i)-0.5;
end

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
                     theta2=-(alphas-alpham2);
                     alpahs_cor=alpham-1.5*theta2;
                elseif itip == 2
                    theta=(alphas-alpham);
                    theta2=(alphas-alpham2);
                    alpahs_cor=alpham+1.5*theta2;
                end
%                 if theta >= 0 && theta < abs(delta)
%                     if itip ==1
%                         alpahs_cor=alpham-theta+delta;
%                     else itip==2
%                         alpahs_cor=alpham+delta+theta;
%                     end
%                 end
                
                if theta > abs(delta)
                    alpahs_cor=-phiw*180/pi+270;
                end
                
                %for hurst method 
                alpahs_x(i,itip)=pi+alpahs_cor*pi/180;  %wave angle to x-axis
                alpahs_N(i,itip)=1.5*pi-(alpahs_cor*pi/180); %wave angle to North
                % refraction coefficient
                kr(itip)=real(sqrt(cos(alpahs_N(i,itip))/cos(phiw)));
                
                % calculate diffraction coeffcient kd from each tip
                if  theta <= 0 && theta >= -90
%                     kd(itip)=0.71-0.0093*theta*(pi/180)+0.000025*(theta*(pi/180))^2;
                    kd(itip)=abs(0.69+0.008*theta*(1));
                    
%                     kd(itip)=0.5*cos(theta*(pi/180)); %hurst
                elseif theta > 0 && theta <= 40
                    kd(itip)=0.71+0.37*sin(theta*(pi/180));
                    %hurst
%                     ls=hypot(Ps2(2,1)-Ps2(2,2),Ps2(1,1)-Ps2(1,2));
%                     xs=hypot(Ps2(2,itip)-y(i),Ps2(1,itip)-x(i));
%                     kd(itip)=0.5*(1-sin(pi*xs/ls));
                elseif theta > 40 && theta <=90
                    kd(itip)=0.83+0.17*sin(theta*(pi/180));
%                     ls=hypot(Ps2(2,1)-Ps2(2,2),Ps2(1,1)-Ps2(1,2));
%                     xs=hypot(Ps2(2,itip)-y(i),Ps2(1,itip)-x(i));
%                     kd(itip)=0.5*(1-sin(pi*xs/ls));
                elseif theta < -90
                    kd(itip)=0;
                else
                    kd(itip)=0;
                end

            end
            
                     %             dp(i)=hbr(i);
            dp(i)=A*(hypot(yp(i)-y(i),xp(i)-x(i))).^(2/3);
            Lp(i) = Lo.*sqrt(tanh(4*pi.^2.*dp(i)./S.tper.^2./S.g));
            
            
            % Shoaling coefficient
            kp(i)=2*pi/Lp(i);
            np(i)=0.5*(1+(2*kp(i)*dp(i)/(sinh(2*kp(i)*dp(i)))));
            
            kI=2*pi/Ltip;
            nI=0.5*(1+(2*kI*dtip/(sinh(2*kI*dtip))));
            
            ks(i_hard,i)=sqrt((nI*Ltip)/(np(i)*Lp(i)));
            
            % combining kd from both tips
                kdt(i_hard,i)=sqrt(kd(1)^2+kd(2)^2+2*kd(1)*kd(2)*0);
            % combining kr from both tips 
%             krt(i_hard,i)=sqrt(kr(1)^2+kr(2)^2+2*kr(1)*kr(2)*0);

          
%        
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             %% combining wave angles from both tips in the shadow zone
            %for hurst method
            hzl_cor=((kd(1)*hstip)^2)*cos(alpahs_x(i,1))+((kd(2)*hstip)^2)*cos(alpahs_x(i,2));
            vl_cor=((kd(1)*hstip)^2)*sin(alpahs_x(i,1))+((kd(2)*hstip)^2)*sin(alpahs_x(i,2));
            
            
            phicor_N=-(atan2(vl_cor,hzl_cor)+0.5*pi)+pi; %from North
            phic(i)=2*pi-atan2(y(i+1)-y(i),x(i+1)-x(i)); %coastline oreintation
            phicor(i)=atan2(sin(phic(i)-phicor_N),cos(phic(i)-phicor_N)); %normal to shoreline

            sphibr(i)=phicor(i);

            tombolo(i)=(mean(ytip)-y(i)) < 1;
            krt(i_hard,i)=real(sqrt(cos(phicor_N)/cos(phiw)));
            
              %% total Coefficients
%                         Kdrs(i_hard,i)=krt(i_hard,i)*kdt(i_hard,i)*ks(i_hard,i);
            Kdrs(i_hard,i)=kdt(i_hard,i)*ks(i_hard,i);
            
        end
    end
end
% plot(xp,yp,'*',x,y,x_mc,y_mc,x_hard,y_hard)
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
    
    if Kdrs(n_hard,i)< 1
        if strcmpi(S.trform,'KAMP')

                    QSkampmass(i)=2.33 * S.rhos/(S.rhos-S.rhow) .* S.tper.^1.5 .* S.tanbeta.^0.75 .* S.d50.^-0.25 .* (hstip*Kdrs(n_hard,i)).^2 .* ( (abs(sin(2*sphibr(i))).^0.6.*sign(sphibr(i)))-(2/S.tanbeta)*cos(sphibr(i)).*(Kdrs(n_hard,ip1)*hstip-Kdrs(n_hard,im1)*hstip)/( hypot(yp(ip1)-yp(im1),(xp(ip1)-xp(im1)))));  %4*S.ds0   
            %         QSkampmass(i)=2.1 * S.rhos/(S.rhos-S.rhow) .* S.tper.^1.5 .* S.tanbeta.^0.75 .* S.d50.^-0.25 .* (hstip*Kdrs(n_hard,i)).^2 .* ( (abs(sin(2*sphibr(i)*phibr_cor(n_hard,i))).^0.6.*sign(sphibr(i)))-(2/S.tanbeta)*cos(sphibr(i)*phibr_cor(n_hard,i)).*(Kdrs(n_hard,ip1)*hstip-Kdrs(n_hard,im1)*hstip)/(2*S.ds0));
            QS(i) = 365*24*60*60*(QSkampmass(i) /(S.rhos-S.rhow)) /(1.0-S.porosity);
        elseif strcmpi(S.trform,'CERC3')
            QS(i)=S.b*(hstip*Kdrs(n_hard,i)).^2.5.*(sin(2*(sphibr(i)))-(2/S.tanbeta)*cos(sphibr(i)).*(Kdrs(n_hard,ip2)-Kdrs(n_hard,im2))/(4*S.ds0));
        end
%     else
%         QS(i)=QS(i);
     QS(tombolo)=0;
    end
end
end


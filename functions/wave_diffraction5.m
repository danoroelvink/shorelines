function [ xS,yS,shadowS ] = find_shadows( x,y,x_mc,y_mc,phiw,hard )

n=length(x)-1;
if hard==1
    crit=0;
else
    crit=0;
end
if n==0
    xS=[];
    yS=[];
    shadowS=[];
else
    f=180/pi;
    xS=.5*(x(1:n)+x(2:n+1));
    yS=.5*(y(1:n)+y(2:n+1));
    len=5*hypot(max(x_mc)-min(x_mc),max(y_mc)-min(y_mc));
    
    for i=1:n
        xw=[xS(i)+1*sin(phiw),xS(i)+len*sin(phiw)];
        yw=[yS(i)+1*cos(phiw),yS(i)+len*cos(phiw)];
        P1=InterX([x_mc;y_mc],[xw;yw]);
        
        shadowS(i)=size(P1,2)>crit;
   shadowS2(i)=size(P1,2)
    end
    
    i=2;
     while i<n
    
       if  shadowS2(i-1) == 0 && shadowS2(i) > 0
           if shadowS2(i) ==2 
                   Hamada=1
           elseif 
                   Hamada=2
               end
               
               xtip(Hamada)=x(i-1);
               ytip(Hamada)=y(i-1);
           while shadowS2(i) > 0
    
              
               %%%%%%%%%%%%%
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

%                                     if yp(i) > mean(ytip)
%                                        yp(i)= mean(ytip);
%                                     end
%% Diffraction Coefficient
            % calculate shadow line angle & theta
            for itip=Hamada
                alpham2=atan2(real((Ps2(2,itip)-ytip(itip))),real(( Ps2(1,itip)-xtip(itip))))*180/pi;
                alpham=atan2(real((Ps1(2,itip)-ytip(itip))),real(( Ps1(1,itip)-xtip(itip))))*180/pi;
                alphas=atan2(yp(i)-(ytip(itip)),(xp(i)-xtip(itip)))*180/pi;
                delta=(alpham-alpham2);
                if itip == 1
                    theta=-(alphas-alpham);
                    alpahs_cor=alphas+delta;
                elseif itip == 2
                    theta=(alphas-alpham);
                    alpahs_cor=alphas+delta;
                end
                if theta >= 0 && theta < abs(delta)
                    if itip ==1
                        alpahs_cor=alpham-theta+delta;
                    else itip==2
                        alpahs_cor=alpham+delta+theta;
                    end
                end
                
                if theta > abs(delta)
                    alpahs_cor=phiw+270;
                end
                
                
                alpahs_x(i,itip)=pi+alpahs_cor*pi/180; %wave angle to x-axis
                alpahs_N(i,itip)=1.5*pi-(alpahs_cor*pi/180); %wave angle to North
                
                % refraction coefficient
                kr(itip)=real(sqrt(cos(alpahs_N(i,itip))/cos(phiw)));
                
                % calculate diffraction coeffcient kd from each tip
                if  theta <= 0 && theta >= -90
                    % kd(itip)=0.71-0.0093*theta*(pi/180)+0.000025*(theta*(pi/180))^2;
                    kd(itip)=abs(0.69+0.008*theta*(1));
                elseif theta > 0 && theta <= 40
                    kd(itip)=0.71+0.37*sin(theta*(pi/180));
                elseif theta > 40 && theta <=90
                    kd(itip)=0.83+0.17*sin(theta*(pi/180));
                elseif theta < -90
                    kd(itip)=0.1;
                else
                    kd(itip)=1;
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
            krt(i_hard,i)=sqrt(kr(1)^2+kr(2)^2+2*kr(1)*kr(2)*0);
            
            
          
            
            %% total Coeffcients
%                         Kdrs(i_hard,i)=krt(i_hard,i)*kdt(i_hard,i)*ks(i_hard,i);
            Kdrs(i_hard,i)=kdt(i_hard,i)*1;
             
               %%%%%%%%%%%%%%%
          i=i+1
           end
       else
           i=i+1
       end
        

    end
    
    
    
    
    
end

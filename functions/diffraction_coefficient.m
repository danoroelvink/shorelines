function [kdt]=diffraction_coefficient(xtip,ytip,x,y,i,Ps1)
for itip=1:2
                alpham=atan2((Ps1(2,itip)-ytip(itip)),( Ps1(1,itip)-xtip(itip)))*180/pi;
                alphas=atan2(y(i)-(ytip(itip)),(x(i)-xtip(itip)))*180/pi;
                if itip == 1
                    theta=-(alphas-alpham);
                elseif itip == 2
                    theta=(alphas-alpham);
                end
                
                % calculate diffraction coeffcient kd from each tip
                if  theta <= 0 && theta >= -90
                   % kd(itip)=0.71-0.0093*theta*(pi/180)+0.000025*(theta*(pi/180))^2;
                   kd(itip)=abs(0.69+0.008*theta*(1));
                elseif theta > 0 && theta <= 40
                    kd(itip)=0.71+0.37*sin(theta*(pi/180));
                elseif theta > 40 && theta <=90
                    kd(itip)=0.83+0.17*sin(theta*(pi/180));
%                 else
%                    kd(itip)=0; 
                end
                
            end
            
%             Lo = S.g*S.tper.^2./2./pi; 
%             L  = Lo.*sqrt(tanh(4*pi.^2.*dtip./S.tper.^2./S.g));
%           
%             shd1=hypot((x_tip(1)-x(i)),(y_tip(1)-y(i)));
%             shd2=hypot((x_tip(2)-x(i)),(y_tip(2)-y(i)));
%             
%             theta_phase=(((shd2-shd1))/L)*360*pi/180;
% %           kdt(i)=sqrt(kd(1)^2+kd(2)^2+2*kd(1)*kd(2)*1);
            kdt(i)=sqrt(kd(1)^2+kd(2)^2+2*kd(1)*kd(2)*0);
          %  kdt(i)=kd(1)*kd(2);%hypot(kd(1),kd(2));
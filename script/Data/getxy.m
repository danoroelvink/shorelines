function [X,Y,TRF] = getxy(xcor,ycor,flipDIM)

    if nargin<3
        flipDIM=[1,1,0];
        TRF.Xoffset=[];
    elseif isstruct(flipDIM)
        TRF=flipDIM;
    else
        TRF.Xoffset=[];
    end
    
    if isempty(TRF.Xoffset)
%         dm1 = diff(xcor,1,1);
%         dn1 = diff(xcor,1,2);
%         dm2 = diff(ycor,1,1);
%         dn2 = diff(ycor,1,2);
%         dm = (dm1.^2+dm2.^2).^0.5;
%         dn = (dn1.^2+dn2.^2).^0.5;
% 
%         X=[zeros(1,size(xcor,2));cumsum(dm,1)];
%         Y=[zeros(size(xcor,1),1),cumsum(dn,2)];

        %% compute X,Y offset, Rotation and FlipDIM
        dx = xcor(end,1)-xcor(1,1);
        dy = ycor(end,1)-ycor(1,1);     
        TRF.Xoffset = xcor(1,1);
        TRF.Yoffset = ycor(1,1); 
        TRF.rota1    = mod(atan2(xcor(end,1)-xcor(1,1),ycor(end,1)-ycor(1,1))*180/pi+180,360)-180;
        TRF.rota2    = mod(atan2(xcor(end,end)-xcor(end,1),ycor(end,end)-ycor(end,1))*180/pi+90+180,360)-180;
        TRF.rota3    = mod(atan2(xcor(1,1)-xcor(1,end),ycor(1,1)-ycor(1,end))*180/pi+270+180,360)-180;
        TRF.rota4    = mod(atan2(xcor(1,end)-xcor(end,end),ycor(1,end)-ycor(end,end))*180/pi+180+180,360)-180;
        TRF.rota = (TRF.rota1+TRF.rota2+TRF.rota3+TRF.rota4)/4+flipDIM(3);
        sa = sin(TRF.rota*pi/180);
        ca = cos(TRF.rota*pi/180);
        X = ca*(xcor-TRF.Xoffset) - sa*(ycor-TRF.Yoffset);
        Y = sa*(xcor-TRF.Xoffset) + ca*(ycor-TRF.Yoffset);

        TRF.flipDIM = flipDIM;
        TRF.Xsize   = 0;
        TRF.Ysize   = 0;
        
        %% flip vertically or horizontally if requested
%         if flipDIM(1)==-1 || flipDIM(3)==-90
%             TRF.Xsize = max(X,[],1);
%             X = repmat(max(X,[],1),[size(xcor,1),1]) - X;
%             %TRF.Ysize = max(Y,[],2);
%         elseif flipDIM(2)==-1
%             TRF.Ysize = max(Y,[],2);
%             Y = repmat(max(Y,[],2),[1,size(xcor,2)]) - Y;
%         end

%         %% flip x and y if requested
%         if abs(flipDIM(3))==90
%             Xtmp = X;
%             Ytmp = Y;
%             X = Ytmp;
%             Y = Xtmp;
%         end
        
    else

        sa = sind(double(-TRF.rota));
        ca = cosd(double(-TRF.rota));
        X = ca.*(double(xcor)+double(TRF.Xoffset)) - sa.*(double(ycor)+double(TRF.Yoffset));
        Y = sa.*(double(xcor)+double(TRF.Xoffset)) + ca.*(double(ycor)+double(TRF.Yoffset));
        
%         sa = sin(-TRF.rota*pi/180);
%         ca = cos(-TRF.rota*pi/180);
%         X = ca*(xcor-TRF.Xoffset) - sa*(ycor-TRF.Yoffset);
%         Y = sa*(xcor-TRF.Xoffset) + ca*(ycor-TRF.Yoffset);

%         if TRF.flipDIM(1)==-1
% %             X = X - mean(TRF.Xsize(:)) * ca;
%             X = X - mean(TRF.Xsize(:));
%             %Y = Y + 2*mean(TRF.Ysize(:));
%         end    
%         if TRF.flipDIM(2)==-1
% %             X = X - mean(TRF.Xsize(:)) * ca;
%             Y = Y + mean(TRF.Ysize(:));
%         end
%          
%         if TRF.flipDIM(3)==90
%             Xtmp = X;
%             Ytmp = Y;
%             X = Ytmp;
%             Y = -Xtmp;
%         end
    end
end

function [structS]=find_struct(x,y,x_hard,y_hard);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
structS=zeros(size(x)-1);
structS=logical(structS);
for i=1:length(x)-1
    P=InterX([[x(i),x(i+1)];[y(i),y(i+1)]],[x_hard;y_hard]);
    structS(i)=size(P,2)>0;
    if 0
        figure(11)
        plot(x,y,'b',[x(i),x(i+1)],[y(i),y(i+1)],'r',x_hard,y_hard,'k',P(1,:),P(2,:),'ro','linewidth',2)
        axis equal
        drawnow
        pause
    end
end
    structS;

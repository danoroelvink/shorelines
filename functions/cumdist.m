function [ s ] = cumdist(x,y)
    s=zeros(size(x));
    for i=2:length(x)
        s(i)=s(i-1)+hypot(x(i)-x(i-1),y(i)-y(i-1));
    end
end


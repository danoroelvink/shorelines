function [xi,yi]=select_multi_polygon(col)
% Select polygon to include in bathy
hold on;
xi = [];yi=[];
n = 0;
% Loop, picking up the points.
finished=0;
while ~finished
    nold=n;
    but = 1;
    while but == 1
        [xs,ys,but] = ginput(1);
        if but==1
            n = n+1;
            xi(n)=xs;
            yi(n)=ys;
            plot(xi,yi,[col '-o']);
        end
    end
    if nold==n || but==113 || but==27                                      % respectively 'q' or 'escape'
        finished=1;
    else
        n=n+1;
        xi(n)=nan;
        yi(n)=nan;
    end
end

% remove single points!
idnan = find(isnan(xi));
xi2=[];
yi2=[];
if length(idnan)>=1
    if idnan(1)>1
        idnan=[0;idnan(:)]';
    end
    if ~isnan(xi(n))
        xi=[xi,nan];
        yi=[yi,nan];
        n=n+1;
        idnan=[idnan,n];
    end
    id = diff(idnan)-1;
    for mm=1:length(id)
        if id(mm)>1
            xi2 = [xi2,xi(idnan(mm)+1:idnan(mm+1))];
            yi2 = [yi2,yi(idnan(mm)+1:idnan(mm+1))];
        end
    end
    xi=xi2;
    yi=yi2;
    n=length(xi);
end

% remove trailing nan
if n>0&&isnan(xi(n))
    xi=xi(1:n-1);
    yi=yi(1:n-1);
end
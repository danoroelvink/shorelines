function [data2]=interpNANs(data,useextrapolation)
% function [data2]=interpNANs(data)
%
% INPUT:
%     data                Array [Nx1] with data
%     useextrapolation    switch to use extrapolation (use 'nearest' or 'linear' for extrapolation)
%                         this does not affect the interpolation which will be 'linear' 
%
% data = [10, 9, nan, 8, nan, 7.5, nan, nan]

if nargin<2
    useextrapolation='nearest';  %<- use extrapolation if it is set to 1
end

rotate=0;
if max(size(data,1),2)<size(data,2)
    rotate=1;
    data=data';
end

if size(data,2)==1
    x = [1:length(data)]';
    y = data;
else
    x = data(:,1);
    y = data(:,2);
end

% make sure that x is increasing
[x,ids]  = sort(x);
y        = y(ids);
eps      = 1e-6;
for tt=1:length(x)-1
    if x(tt)==x(tt+1)
        x(tt+1)=x(tt+1)+eps;
    end
end

%% remove nans
id = find(~isnan(y));
x2 = x(id);
y2 = y(id);

if length(id)==length(x)
     ynew = y2;
else
    if length(x2)>1
        if strcmpi(useextrapolation,'linear')
            ynew = interp1(x2,y2,x,'linear','extrap');
        else
            y3 = interp1(x2,y2,x,'linear');
            id = find(~isnan(y3));
            x2 = x(id);
            y2 = y3(id);
            ynew = interp1(x2,y2,x,'nearest','extrap');
        end
    elseif length(x2)==1
        ynew=repmat(y2,size(x));
    else
        %fprintf('Need more data points than 1\n')
        ynew=nan(size(x));
    end
end

if size(data,2)==2
    data2 = [x,ynew];
else
    data2 = [ynew];
end

if rotate==1
    data2 = data2';
end



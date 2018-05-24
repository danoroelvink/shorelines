function [data2]=interpNANsDIR(data)
% function [data2]=interpNANsDIR(data)
%
% INPUT:
%     data                Array [Nx1] with data
%     useextrapolation    switch to use extrapolation (use 'nearest' or 'linear' for extrapolation)
%                         this does not affect the interpolation which will be 'linear' 
%
% [data2]=interpNANsDIR([nan;nan;30;nan;330;300;nan;320;300;nan;nan;nan;20;nan;nan;350;30;nan;285]);

rotate=0;
if max(size(data,1),2)<size(data,2)
    rotate=1;
    data=data';
end

data1 = data;
dir0 = data(find(~isnan(data),1));
for ii=1:length(data)
    if ~isnan(data(ii));
        dir = data(ii);
        diroffset = dir0-mod(dir0,360);
        dir0B = mod(dir0+180,360);
        if dir<dir0B && abs(dir+diroffset-dir0)>180
            diroffset = diroffset+360;
        elseif dir>dir0B && abs(dir+diroffset-dir0)>180
            diroffset = diroffset-360;
        end
        data1(ii) = dir+diroffset;
        dir0 = data1(ii);
    end
end

data2 = mod(interpNANs(data1,'nearest'),360);

if rotate==1
    data2 = data2';
end



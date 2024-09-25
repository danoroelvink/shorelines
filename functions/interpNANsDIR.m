function [data2]=interpNANsDIR(data)
% function [data2]=interpNANsDIR(data)
%
% This function interpolates NaNs in vectors with directional information.
% For extrapolation the nearest value is used. 
%
% INPUT: 
%     data           : Array [Nx1] with data
%
% [data2]=interpNANsDIR([nan;nan;30;nan;330;300;nan;320;300;nan;nan;nan;20;nan;nan;350;30;nan;285]);
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2020 IHE Delft & Deltares
%
%       Dano Roelvink
%       d.roelvink@un-ihe.org
%       Westvest 7
%       2611AX Delft
%
%       Bas Huisman
%       bas.huisman@deltares.nl
%       Boussinesqweg 1
%       2629HV Delft
%
%   This library is free software: you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation, either
%   version 2.1 of the License, or (at your option) any later version.
%
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses>
%   --------------------------------------------------------------------

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
end

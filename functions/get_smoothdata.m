function [VAR]=get_smoothdata(VAR,smoothtype,smoothsteps,movwindow)
% function [VAR]=get_smoothdata(VAR,smoothtype,smoothsteps,movwindow)
% 
% INPUT: 
%     VAR         : column with the variable [1xN] 
%     smoothtype  : e.g. 'angle'/'vector' or 'scalar' 
%     smoothsteps : number of times the smaoothing is applied (with 'angle' also a fraction between 0 and 1 can be used to apply an alfa) 
%     movwindow   : (optional) using a moving average instead of 'central averaging' 
% 
% OUTPUT:
%     VAR         : smoothed variable [1xN] 
% 
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2021 IHE Delft & Deltares
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
    
    if nargin<3
        smoothsteps=1;
        movwindow=2;
    end

    if isempty(VAR)
       smoothsteps=0;
    end

    if size(VAR,1)>1 && size(VAR,2)==1
        VAR=VAR(:)';
    end

    for kk=1:size(VAR,1)
        VAR2=VAR(kk,:);
        for nn=1:max(smoothsteps,1)
            if strcmpi(smoothtype,'avgmean') 
                VAR2=mean(VAR2);        
            elseif ((strcmpi(smoothtype,'angle') ||  strcmpi(smoothtype,'vector')) && nargin<4) ...
                   || strcmpi(smoothtype,'anglemean') && nn<max(smoothsteps,1)
                sVAR=sind(VAR2(:))';
                cVAR=cosd(VAR2(:))';
                sVARavg=[sVAR(1),(sVAR(1:end-1)+sVAR(2:end))/2,sVAR(end)];
                cVARavg=[cVAR(1),(cVAR(1:end-1)+cVAR(2:end))/2,cVAR(end)];
                sVAR2=(sVARavg(1:end-1)+sVARavg(2:end))/2;
                cVAR2=(cVARavg(1:end-1)+cVARavg(2:end))/2;
                VAR2=atan2d(sVAR2,cVAR2);
                VAR2=mod(VAR2,360);   
                if smoothsteps<1
                    alfa=smoothsteps;
                    VAR2=atan2d((1-alfa)*sVAR+alfa*sVAR2,(1-alfa)*cVAR+alfa*cVAR2);
                    VAR2=mod(VAR2,360);
                end 
            elseif (strcmpi(smoothtype,'angle') || strcmpi(smoothtype,'vector')) && nargin==4
                cVAR=smoothdata(cosd(VAR2),'movmean',movwindow);
                sVAR=smoothdata(sind(VAR2),'movmean',movwindow);
                VAR2=atan2d(sVAR,cVAR);
                VAR2=mod(VAR2,360);
            elseif strcmpi(smoothtype,'angleavgmean')
                sVAR=mean(sind(VAR2));
                cVAR=mean(cosd(VAR2));
                VAR2=atan2d(sVAR,cVAR);
                VAR2=mod(VAR2,360);  
            elseif nargin==4
                VAR2=smoothdata(VAR2,'movmean',movwindow);
            elseif strcmpi(smoothtype,'anglemean') && nargin<4
                sVAR=sind(VAR2(:))';
                cVAR=cosd(VAR2(:))';
                sVARavg=(sVAR(1:end-1)+sVAR(2:end))/2;
                cVARavg=(cVAR(1:end-1)+cVAR(2:end))/2;
                VAR2=atan2d(sVARavg,cVARavg);
                VAR2=mod(VAR2,360);   
            else
                VAR2=VAR2(:)';
                VARavg=[VAR2(1),(VAR2(1:end-1)+VAR2(2:end))/2,VAR2(end)];
                VAR2=(VARavg(1:end-1)+VARavg(2:end))/2;  
            end
        end
        if ~(strcmpi(smoothtype,'anglemean') && nn==max(smoothsteps,1))
            VAR(kk,:)=VAR2;
        else
            VAR(kk,1:end-1)=VAR2;
            if kk==size(VAR,1)
                VAR=VAR(:,1:end-1);
            end
        end
    end
end

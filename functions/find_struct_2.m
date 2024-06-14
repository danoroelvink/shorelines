function [structS]=find_struct_2(x,y,x_hard,y_hard);
% function [structS]=find_struct_2(x,y,x_hard,y_hard);
%
% UNTITLED Summary of this function goes here
% Detailed explanation goes here
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
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses
%   --------------------------------------------------------------------

    
    if (length(x)~=length(y))
        error('find_struct:: length of x is not equal to length of y.');
    end
    if (length(x_hard)~=length(y_hard))
        error('find_struct:: length of x_hard is not equal to length of y_hard.');
    end
    structS=zeros(1,length(x)-1);
    structS=logical(structS);
    for i=1:length(x)-1  
        [xx,yy]=get_intersections([x(i),x(i+1)],[y(i),y(i+1)],x_hard,y_hard);
        if~isempty(xx)
            for j = 1 :length(xx)
                [~, indx] = min(hypot([x(i),x(i+1)]-xx(j),[y(i),y(i+1)]-yy(j)));
                if indx==1
                    structS(i)=1;
                elseif indx==2
                    structS(i+1)=1;
                end
                if 0
                    figure(11)
                    plot(x,y,'b',[x(i),x(i+1)],[y(i),y(i+1)],'r',x_hard,y_hard,'k',xx,yy,'ro','linewidth',2)
                    axis equal
                    drawnow
                    %%pause
                end
            end
        end
    end
end    
    
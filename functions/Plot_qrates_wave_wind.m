function [qwave,qwind,innt,time,step]=plot_qrates_wave_wind(S,x_dune,y_dune,x_trans,y_trans,n_trans,qss,qww,innt,tnow,qwave,qwind,time,step)
% function [qwave,qwind,innt,time,step]=plot_qrates_wave_wind(S,x_dune,y_dune,x_trans,y_trans,n_trans,qss,qww,innt,tnow,qwave,qwind,time,step)
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

    
    for i =1:n_trans
        xs_trans=[x_trans(2*i-1) x_trans((2*i-1)+1)];
        ys_trans=[y_trans(2*i-1) y_trans((2*i-1)+1)];
        P1=InterX([xs_trans;ys_trans],[x_dune;y_dune]);
        [P1,P2]=get_intersections(xs_trans,ys_trans,x_dune,y_dune);
        [~, indx] = min(hypot(x_dune-P1,y_dune-P2));
        qwave(i,innt)=qss(indx);
        qwind(i,innt)=qww(indx);         
    end
    time(innt)=tnow;
    step(innt)=S.dt;
    innt=innt+1;
end

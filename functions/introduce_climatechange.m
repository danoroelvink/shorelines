function [CC]=introduce_climatechange(CC,TIME,i_mc)
% function [CC]=introduce_climatechange(CC,TIME,i_mc)
%
% Derive climate change related corrections at considered time instance.
%
% INPUT: 
%    CC
%        .timenum  : time in datenum
%        .SLR      : rate of sea level rise [m/yr]
%        .HS       : change in wave height since start simulation (m)
%        .DIR      : change in wave direction since start of simulation (deg)
%
%    TIME          : current moment in time 'tnow' is used (in datnum format [days since 0000-00-01 00:00])
%        .timenum0 : time at model start t0 [in days in datenum format]
%        .tnow     : current time [in days in datenum format]
%        .tprev    : time at previous coastal timestep [in days in datenum format]
%
% OUTPUT:
%    CC
%        .SLRo     : rate of sea level rise at considered time instance [m/yr]
%        .HScor    : relative increase of the significant wave height at considered time instance w.r.t. t0 [-]
%        .PHIcor   : relative change in wave direction at considered time instance w.r.t. t0 [delta °]
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2022 IHE Delft & Deltares
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

    if i_mc==1
        if isnan(CC.timenum); return; end

        %% Constant rate per year
        if isempty(CC.timenum)
            % elapsed time
            dt=(TIME.tnow-TIME.timenum0)/365.;
            CC.SLRo   = CC.SLR;              % rate of sea level rise at considered time instance (m/yr)
            CC.HScor  = 1 + CC.HS * dt;      % relative increase of Hs at considered moment w.r.t. t0 [-]
            CC.PHIcor = CC.DIR * dt;         % rotation of the wave direction at considered moment w.r.t. t0 [delta °]

        %% Time series    
        else
            dt      = (TIME.tnow-TIME.tprev)/365.;

            % some checks of the time -> converting strings of years to datenum
            if mean(CC.timenum)>1e7 % time as a string of yearmonthday's
                CC.timenum=datenum(num2str(CC.timenum,'%1.0f'),'yyyymmdd');
            elseif mean(CC.timenum)<1e4 % time as a string of years
                CC.timenum=datenum(num2str(CC.timenum,'%1.0f'),'yyyy');
            end

            if length(CC.timenum(:))==length(CC.SLR(:))
                % within prescribed interval
                if dt>0.0 && TIME.tprev>=min(CC.timenum(:)) && TIME.tnow<=max(CC.timenum(:))
                    slr     = interp1(CC.timenum(:),CC.SLR(:),TIME.tnow);   % SLR amount since simulation start in [m]
                    slrprev = interp1(CC.timenum(:),CC.SLR(:),max(TIME.tprev,min(CC.timenum)));
                    CC.SLRo = (slr-slrprev)/dt;   % instantanious rate for coastline_change

                % before time-series start take initial slr-rate,
                % after time-series take last slr-rate
                else
                    slrrate = diff(CC.SLR(:))./diff(CC.timenum(:));
                    if TIME.tnow<=min(CC.timenum(:)) || TIME.tprev<=min(CC.timenum(:)) || dt==0
                        CC.SLRo = slrrate(1);   % instantanious rate for coastline_change
                    else
                        CC.SLRo = slrrate(end);   % instantanious rate for coastline_change
                    end
                end

            else    
                CC.SLRo=0.0;
            end

            if length(CC.timenum(:))==length(CC.HS(:))
                % CC.HS is the relative increase in wave height w.r.t. HS at t0 [increase as fraction of initial wave height] (e.g. 0.04 means an increase of 4% in wave height w.r.t. t0)
                hsfactor     = interp1(CC.timenum(:),CC.HS(:),TIME.tnow);  
                if isnan(hsfactor)
                    error('Could not determine slr related wave height change factor. Check your input time series');
                end
                CC.HScor = 1+hsfactor;
            end

            if length(CC.timenum(:))==length(CC.DIR(:))
                sinrate      = interp1(CC.timenum(:),sind(CC.DIR(:)),TIME.tnow);
                cosrate      = interp1(CC.timenum(:),cosd(CC.DIR(:)),TIME.tnow);
                dirchange    = atan2d(sinrate,cosrate);
                if isnan(dirchange)
                    error('Could not determine the absolute change in wave direction. Check your input time series');
                end
                CC.PHIcor = dirchange;
            end
        end
    end
end  

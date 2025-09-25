function [wberm,qs_avg,qss_avg,ql_avg,R]=dune_erosion(runupform,runupfactor,wberm_old,dfelev0,dfelev,dcelev,h0,sqrHsL0,SWL,TPc,cs,cstill,xtill,perctill,A,dt,dtnum)
% function [wberm,qs_avg,qss_avg,ql_avg,R]=dune_erosion(runupform,runupfactor,wberm_old,dfelev,dcelev,h0,sqrHsL0,SWL,TPc,cs,cstill,xtill,perctill,A,dt,dtnum)
% 
% The rate of erosion of the dune is computed in this routine with the
% formulation of Larson et al. on the basis of the wave conditions 
% and runup height/still water level. 
% The change in the berm width is computed as well as the erosion from
% the dunes (qs+ql) and the sediment that becomes available for the
% beach to grow. 
% 
% INPUT:
%    runupform       : runup formulation applied
%    runupfactor     : tuning factor for runup, considering inaccuracies in runup formulations
%    wberm_old       : berm width (distance MSL to dune foot) [m]
%    dfelev0         : dune foot elevation, original used for assessing the active height [m]
%    dfelev          : dune foot elevation, adjusted for beach width [m] (i.e. with a reduced dunefoot when the beach is too narrow, meaning steeper than the maxslope)
%    dcelev          : dune crest elevation [m]
%    h0              : active height of the profile [m]
%    sqrHsL0         : square of the wave height and wave length
%    SWL             : still water level [m]
%    TPc             : wave period [s]
%    cs              : erosion coefficient of the dunes during storms, which scales the rate of erosion
%    cstill          : erosion coefficient of the dunes during storms for dunes with consolidated till layers, which scales the rate of erosion
%    xtill           : width of sandy dune in front of cohesive dune
%    perctill        : percentage of till (0 to 100) of the cohesive dune, with the sand percentage being 100-perctill
%    A               : overwash parameter A
%    dt              : dune timestep [year]
%    dtnum           : number of wind/wave/swl conditions within a single coastline timestep [-]
% 
% OUTPUT:
%    wberm           : berm width (distance MSL to dune foot) [m]
%    qs_avg          : dune erosion volume change that increases the beach width [m3/m/yr]
%    qss_avg         : dune erosion volume change that increases the beach width, sandy part only [m3/m/yr]
%    ql_avg          : dune erosion volume change that does not increase the beach width [m3/m/yr]
%    R               : runup level for dune erosion [m]
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
   
   wberm=wberm_old;
   spyr=365*24*3600;
   nt=round(dt/dtnum,0);
   nt=max(nt,1);
   qd=zeros(size(wberm));
   qs=zeros(size(wberm));
   qss=zeros(size(wberm));
   ql=zeros(size(wberm));
   alfa=zeros(size(wberm));
   qs_avg=zeros(size(wberm));
   qss_avg=zeros(size(wberm));
   ql_avg=zeros(size(wberm));
   dcelev=max(dcelev,dfelev0+0.2);
   DH=dcelev-dfelev0;
   wbermmin=5;
   for step=1:nt
        slope=max(min(dfelev./(wberm-wbermmin),0.5),0.001);
        if ~isempty(findstr(lower(runupform),'sto'))
            % Run-up Height [m] Stockdon et al. (2006)
            R = runupfactor*1.1*sqrHsL0.*(0.35*slope + sqrt(0.563*(slope.^2) + 0.004)/2);
        elseif ~isempty(findstr(lower(runupform),'gho'))
            % Run-up Height [m] thesis M. Ghonim
            R = runupfactor*slope.*sqrHsL0;
        else
            % Run-up height [m] Larson et al. (2016) = default
            %CF = 0.03;
            %X = 2 .* wberm .* (1 - SWL ./ dfelev);
            R = runupfactor*0.158*sqrHsL0;
            %RPRIM = R .* exp(-2. .* CF .* X) + (dfelev - SWL) .* (1. - exp(-2. .* CF .* X));
            %if RPRIM + SWL > DFOOT
            %    R = RPRIM;
            %end
        end
        
        ero=SWL+R>dfelev;
        wldf=SWL+R-dfelev;
        if sum(ero)==0
            break
        end
        
        % maximum avaiable volume to be eroded from sandy part of the dune
        % erosion beyond xtill it will have cohesive behaviour
        % by default assume unlimited width of the sandy dunes (1e6) when till is not defined. 
        qdmax=5.0*DH/dtnum; 
        if ~isempty(xtill)
        qdmax=xtill.*DH/dtnum;
        end
        
        % Dune erosion flux according to Larson (??) for the sandy part
        % with qd in m3/m/yr
        qd(ero)=spyr * 4.*cs(ero).*(wldf(ero)).* min(wldf(ero),DH(ero))./ TPc(ero);
        qd(ero)=min(qd(ero),qdmax(ero));  % To limit erosion/overwash per dtnum
        qds=qd;
        
        % Dune erosion flux that is cohesive
        tillfrac=max(1-qdmax(ero)./qd(ero),0); % fraction of erosion that is not anymore available as sand in the dune, but should come from the cohesive till layer behind the sand
        qd2=spyr * tillfrac.*4.*cstill(ero).*(wldf(ero)).* min(wldf(ero),DH(ero))./ TPc(ero);
        qd(ero)=qd(ero)+qd2;
        qds(ero)=qd(ero)+min(1,max(1-perctill(ero)/100,0)).*qd2;
        
        % differentiate between the erosion components qs and ql
        alfa(ero)=max(((wldf(ero))./DH(ero)-1),0)/A;
        qs(ero)=qd(ero)./(1+alfa(ero));
        qss(ero)=qds(ero)./(1+alfa(ero));
        ql(ero)=qd(ero).*alfa(ero)./(1+alfa(ero));
        qs(~ero)=0;
        qd(~ero)=0;
        qds(~ero)=0;
        qss(~ero)=0;
        ql(~ero)=0;

        % change in berm width (wberm) as a result of changes in coastline position (dn) and dunefront position (dndune)
        % The qds is the sand that is eroded from the dune, which is the same as the total erosion for sandy dunes. 
        % In case of cohesive till layers the 'perctill' will determine how much sediment (of the erodeded dune material) becomes available to the coast/beach. 
        dn=qds./h0*dtnum;
        dndune=(qs+ql)./DH*dtnum;
        wberm=max(wberm+dn+dndune,0.1);
        
        % average erosion/accretio rate
        qs_avg=qs_avg+qs;
        qss_avg=qss_avg+qss;
        ql_avg=ql_avg+ql;
    end
    qs_avg=qs_avg/nt;
    qss_avg=qss_avg/nt;
    ql_avg=ql_avg/nt;
end


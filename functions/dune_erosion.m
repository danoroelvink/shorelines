function [Wberm,qs_avg,qss_avg,ql_avg,R,xtill]=dune_erosion(runupform,runupfactor,Wberm_old,Dfelev,Dcelev,h0,sqrHsL0,SWL,TPc,Cs,Cstill,xtill,perctill,A,dt,dtnum)
% function [Wberm,qs_avg,qss_avg,ql_avg,R,xtill]=dune_erosion(runupform,runupfactor,Wberm_old,Dfelev,Dcelev,h0,sqrHsL0,SWL,TPc,Cs,Cstill,xtill,perctill,A,dt,dtnum)
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
   
   Wberm=Wberm_old;
   spyr=365*24*3600;
   nt=round(dt/dtnum,0);
   nt=max(nt,1);
   dtact=dt/nt;
   qd=zeros(size(Wberm));
   qs=zeros(size(Wberm));
   qss=zeros(size(Wberm));
   ql=zeros(size(Wberm));
   alfa=zeros(size(Wberm));
   qs_avg=zeros(size(Wberm));
   qss_avg=zeros(size(Wberm));
   ql_avg=zeros(size(Wberm));
   Dcelev=max(Dcelev,Dfelev+0.2);
   DH=Dcelev-Dfelev;
   for step=1:nt
        slope=min(Dfelev./Wberm,0.5);
        if ~isempty(findstr(lower(runupform),'sto'))
            % Run-up Height [m] Stockdon et al. (2006)
            R = runupfactor*1.1*sqrHsL0.*(0.35*slope + sqrt(0.563*(slope.^2) + 0.004)/2);
        elseif ~isempty(findstr(lower(runupform),'gho'))
            % Run-up Height [m] thesis M. Ghonim
            R = runupfactor*slope.*sqrHsL0;
        else
            % Run-up height [m] Larson et al. (2016) = default
            R = runupfactor*0.158*sqrHsL0;
        end
        ero=R+SWL>Dfelev;
        if sum(ero)==0
            break
        end
        
        % maximum avaiable volume to be eroded from sandy part of the dune
        % erosion beyond xtill it will have cohesive behaviour
        qdmax=xtill.*h0.*dt;

        % Dune erosion flux according to Larson (??) for the sandy part
        % with qd in m3/m/yr
        qd(ero)=spyr * 4.*Cs(ero).*(SWL(ero)+R(ero)-Dfelev(ero)).* ...
         min(SWL(ero)+R(ero)-Dfelev(ero),DH(ero))./ TPc(ero);
        qd(ero)=min(qd(ero),5*DH(ero)/dtnum);  % To limit erosion/overwash per dtnum
        qds=qd;

        % Dune erosion flux that is cohesive
        qd1=min(qd(ero),qdmax(ero));
        tillfrac=max(1-qdmax(ero)./qd(ero),0); % fraction of erosion that is not anymore available as sand in the dune, but should come from the cohesive till layer behind the sand
        qd2=spyr * tillfrac.*4.*Cstill(ero).*(SWL(ero)+R(ero)-Dfelev(ero)).* ...
          min(SWL(ero)+R(ero)-Dfelev(ero),DH(ero))./ TPc(ero);
        qd(ero)=qd1+qd2;
        qds(ero)=qd1+min(1,max(1-perctill(ero)/100,0)).*qd2;

        % differentiate between the erosion components qs and ql
        alfa(ero)=max(((SWL(ero)+R(ero)-Dfelev(ero))./DH(ero)-1),0)/A;
        qs(ero)=qd(ero)./(1+alfa(ero));
        qss(ero)=qds(ero)./(1+alfa(ero));
        ql(ero)=qd(ero).*alfa(ero)./(1+alfa(ero));
        qs(~ero)=0;
        qd(~ero)=0;
        qds(~ero)=0;
        qss(~ero)=0;
        ql(~ero)=0;

        % change in berm width (Wberm) as a result of changes in coastline position (dn) and dunefront position (dndune)
        % The qds is the sand that is eroded from the dune, which is the same as the total erosion for sandy dunes. 
        % In case of cohesive till layers the 'perctill' will determine how much sediment (of the erodeded dune material) becomes available to the coast/beach. 
        dn=qds./h0*dtnum;
        dndune=(qs+ql)./DH*dtnum;
        Wberm=Wberm+dn+dndune;

        % average erosion/accretio rate
        qs_avg=qs_avg+qs;
        qss_avg=qss_avg+qss;
        ql_avg=ql_avg+ql;
    end
    qs_avg=qs_avg/nt;
    qss_avg=qss_avg/nt;
    ql_avg=ql_avg/nt;
end


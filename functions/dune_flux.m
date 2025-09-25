function [DUNE,COAST]=dune_flux(COAST,DUNE,WIND,RUNUP,TRANSP,TIME,CC)
% function [DUNE,COAST]=dune_flux(COAST,DUNE,WIND,RUNUP,TRANSP,TIME,CC)
% 
% This routine computes the dune erosion and dune growth.
% 
%  COAST
%        .h0             : active height of the profile [m]
%        .PHIcxy         : shore-normal orientation at coastline points [°N]
%        .wberm          : Berm width (distance MSL to dune foot) [m]
%        .dfelev         : Dune foot elevation [m]
%        .dcelev         : Dune crest elevation [m]
%        .cs             : Erosion coefficient of the dunes during storms, which scales the rate of erosion
%        .cstill         : Erosion coefficient of the dunes during storms for dunes with consolidated till layers, which scales the rate of erosion
%        .xtill          : Width of sandy dune in front of cohesive dune
%        .perctill       : Percentage of till (0 to 100) of the cohesive dune, with the sand percentage being 100-perctill
%  WIND
%        .uz             : wind velocity [m/s]
%        .phiwnd         : wind direction [°N]
%        .rhoa           : Density of air [kg/m3]
%  RUNUP
%        .swl            : surge water-level [m w.r.t. MSL]
%        .Hs             : significant wave height [m]
%        .Tp             : wave period [s]
%  DUNE
%        .dt             : dune timestep [year]
%        .A              : Overwash parameter A
%        .duneaw         : Coefficient (Bagnold, 1937)
%        .rhos           : Density of the sediment [kg/m3]
%        .z              : wind measurement vertical height [m]
%        .g              : Gravitational acceleration [m/s2]
%        .d50            : Median grain diameter [m]
%        .d50r           : Median reference grain size [m]
%        .k              : Von Karman's coefficient
%        .kw             : Empirical coefficient (Sherman et al. 2013)
%        .porosity       : porosity
%        .segmaw         : Empirical factor used for scaling impact of the fetch length
%        .maxslope       : The maximum slope angle of the dunes (1:slope). The dunefoot height is lowered if the beach gets too steep (preserving the max slope).
%        .runupform      : runup formulation applied
%        .runupfactor    : tuning factor for runup, considering inaccuracies in runup formulations
%  TRANSP
%        .shadow         : index of cells with shadowing of other parts of the coastlines (xy-points)
%        .shadow_h       : index of cells with shadowing of hard structures (xy-points)
%        .shadowS_hD     : index of cells with shadowing of hard structures (QS-points)
%  TIME
%        .dt             : time step [year]
% 
% OUTPUT: 
%  COAST
%        .wberm          : width of the berm/beach [m]
%        .qs             : dune erosion volume change that increases the beach width [m3/m/yr]
%        .qss            : dune erosion volume change that increases the beach width, sandy part only [m3/m/yr]
%        .ql             : dune erosion volume change that does not increase the beach width [m3/m/yr]
%        .qw             : wind transport from the beach to the dune [m3/m/yr]
%        .R              : runup level for dune erosion [m]
%        .SWL            : still water level for dune erosion [m]
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
    
    %% Built-in dune model of ShorelineS
    if DUNE.used
        dt       = TIME.dt;
        if ~isempty(DUNE.dt)
        dtdune   = DUNE.dt;
        else
        dtdune   = TIME.dt;
        end
        dtnum    = min(dtdune,1/365/48);
        phic     = COAST.PHIcxy;
        h0       = COAST.h0;
        n        = length(phic);
        phiwnd   = WIND.phiwnd;
        uwind    = WIND.uz;
        rhoa     = WIND.rhoa;
        SWL      = RUNUP.swl;
        HS       = RUNUP.Hs;
        TP       = RUNUP.Tp;
        dPHI     = repmat(COAST.PHIcxy,[size(RUNUP.Dir,1),1])-RUNUP.Dir;
        xtill    = COAST.xtill;
        perctill = COAST.perctill;
        cs       = COAST.cs(:)';
        cstill   = COAST.cstill(:)';
        A        = DUNE.A;
        duneaw   = DUNE.duneaw;
        z        = DUNE.z;
        rhos     = DUNE.rhos;
        g        = DUNE.g;
        d50      = DUNE.d50;
        d50r     = DUNE.d50r;
        k        = DUNE.k;
        kw       = DUNE.kw;
        por      = DUNE.porosity;
        segmaw   = DUNE.segmaw;
        maxslope = DUNE.maxslope;
        runupform= DUNE.runupform;
        runupfactor=DUNE.runupfactor;
        
        %% Berm width
        wberm    = COAST.wberm;
        dfelev0  = COAST.dfelev;
        dcelev   = COAST.dcelev;
        
        % Extend single values over the whole grid
        %uwind=repmat(uwind,[1,n/size(uwind,2)]);
        %phiwnd=repmat(phiwnd,[1,n/size(phiwnd,2)]);
        
        %% Initialize fluxes
        COAST.qs  = zeros(1,size(COAST.x,2));
        COAST.qss = zeros(1,size(COAST.x,2));
        COAST.qw  = zeros(1,size(COAST.x,2));
        COAST.ql  = zeros(1,size(COAST.x,2));
        COAST.R   = zeros(1,size(COAST.x,2));
        COAST.SWL = zeros(1,size(COAST.x,2));

        % make sure to use sub-timesteps when multiple measurements of HS, SWL, TP or uwind are taking place in a single timestep 
        % (i.e. at more than just one timestep, but also inbetween moments)
        dtfractions=max([size(HS,1),size(TP,1),size(dPHI,1),size(SWL,1),size(uwind,1),size(phiwnd,1)]);
        
        % Loop over the sub-conditions of the time step (i.e. at fractions of coastline timestep)
        for tt=1:dtfractions
            HSi=HS(min(tt,size(HS,1)),:);
            TPi=TP(min(tt,size(TP,1)),:);
            dPHIi=dPHI(min(tt,size(dPHI,1)),:);
            SWLi=SWL(min(tt,size(SWL,1)),:);
            uwindi=uwind(min(tt,size(uwind,1)),:);
            phiwndi=phiwnd(min(tt,size(phiwnd,1)),:);
            
            % time step
            dtdune=dt/dtfractions;
            
            % Slope of the berm
            % make sure to limit the beach slope, which effectively means the dunefoot will get lower for narrow beaches, and the odds of erosion will increase.
            slope0=dfelev0./wberm;
            dfelev=dfelev0.*min(maxslope./slope0,1);
            slope=dfelev./wberm;
            
            % Offshore wave conditions (or at toe of dynamic profile(?)
            HSo1 = HSi.*sqrt(max(cosd(dPHIi),eps));
            L0 = g/(2*pi)*TPi.^2;
            sqrHsL0=sqrt(HSo1.*L0);
            
            %% Compute dune erosion flux by waves
            %SWLi=get_smoothdata(SWLi,'',ceil(250/COAST.ds0));
            %TPi=get_smoothdata(TPi,'',ceil(250/COAST.ds0));
            %sqrHsL0=get_smoothdata(sqrHsL0,'',ceil(250/COAST.ds0));
            [~,qs,qss,ql,R]=dune_erosion(runupform,runupfactor,wberm,dfelev0,dfelev,dcelev,h0,sqrHsL0,SWLi,TPi,cs,cstill,xtill,perctill,A,dtdune,dtnum);
            
            %% Aeolian transport rate from beach to dune
            [qw,uc,us,mwe,qwe,Bdry,FL]=dune_growth(uwindi,duneaw,rhos,rhoa,por,g,d50,d50r,z,k,kw,R,SWLi,dfelev,wberm,phic,phiwndi,segmaw);
            
            % account for probability of the condition
            Prob=1/dtfractions; % in case of eqaul likelihood (e.g. time-series)
            if isfield(WIND,'Prob') && isfield(WIND,'iwc')
                Prob=WIND.Prob(WIND.iwc(tt))/sum(WIND.Prob(WIND.iwc)); % in case of varying probability (e.g. wave climate)
            end
            
            % set to 0 for shadowed sections of coast
            tt2=min(tt,size(TRANSP.shadow,1));
            qs(TRANSP.shadow(tt2,:))      = 0;      % Set dune transport to zero due to dune-dune shadowing
            qs(TRANSP.shadow_h(tt2,:))    = 0;      % Set dune transport to zero due to being sheltered behind structures
            qs(TRANSP.shadowS_hD(tt2,:))  = 0;      % Set dune transport to zero due to being sheltered behind revetment
            qss(TRANSP.shadow(tt2,:))     = 0;      % Set dune transport to zero due to dune-dune shadowing
            qss(TRANSP.shadow_h(tt2,:))   = 0;      % Set dune transport to zero due to being sheltered behind structures
            qss(TRANSP.shadowS_hD(tt2,:)) = 0;      % Set dune transport to zero due to being sheltered behind revetment
            ql(TRANSP.shadow(tt2,:))      = 0;      % Set dune transport to zero due to dune-dune shadowing
            ql(TRANSP.shadow_h(tt2,:))    = 0;      % Set dune transport to zero due to being sheltered behind structures
            ql(TRANSP.shadowS_hD(tt2,:))  = 0;      % Set dune transport to zero due to being sheltered behind revetment   
            
            % export data
            COAST.qw=COAST.qw+qw*Prob;         % compute average qw, i.e. component due to wind transport
            COAST.qs=COAST.qs+qs*Prob;         % compute average qs, i.e. dune erosion that is beneficial for the beach 
            COAST.qss=COAST.qss+qss*Prob;      % compute average qss, i.e. dune erosion that is beneficial for the beach (i.e. the sandy part, discarding the till part)
            COAST.ql=COAST.ql+ql*Prob;         % compute average ql, i.e. dune erosion that is not beneficial for the beach 
            COAST.R=max(COAST.R,R);            % use maximum instead of average for the runup
            COAST.SWL=max(COAST.SWL,SWLi);     % use maximum instead of average for the still water level
        end
    end
end

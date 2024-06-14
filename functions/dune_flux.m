function [DUNE,COAST]=dune_flux(COAST,DUNE,WIND,RUNUP,TRANSP,TIME)
% function [DUNE,COAST]=dune_flux(COAST,DUNE,WIND,RUNUP,TRANSP,TIME)
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

    if DUNE.used
        dt       = TIME.dt;
        dtdune   = DUNE.dt;
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
        Cs       = COAST.Cs(:)';
        Cstill   = COAST.Cstill(:)';
        A        = DUNE.A;
        duneAw   = DUNE.duneAw;
        z        = DUNE.z;
        rhos     = DUNE.rhos;
        g        = DUNE.g;
        d50      = DUNE.d50;
        d50r     = DUNE.d50r;
        k        = DUNE.k;
        Kw       = DUNE.Kw;
        por      = DUNE.porosity;
        segmaw   = DUNE.segmaw;
        maxslope = DUNE.maxslope;
        spyr     = 3600*24*365;
        dtnum    = 1/365/48;
        runupform= DUNE.runupform;
        runupfactor=DUNE.runupfactor;
        
        %% Berm width
        Wberm    = COAST.Wberm;
        Dfelev0  = COAST.Dfelev;
        Dcelev   = COAST.Dcelev;
        
        % Extend single values over the whole grid
        uwind=repmat(uwind,[1,n/size(uwind,2)]);
        phiwnd=repmat(phiwnd,[1,n/size(phiwnd,2)]);
        
        %% Initialize fluxes
        COAST.qs  = zeros(size(COAST.x));
        COAST.qss = zeros(size(COAST.x));
        COAST.qw  = zeros(size(COAST.x));
        COAST.ql  = zeros(size(COAST.x));
        COAST.R   = zeros(size(COAST.x));
        COAST.SWL = zeros(size(COAST.x));

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
            
            %% Compute dune erosion flux by waves
        
            % Slope of the berm
            % make sure to limit the beach slope, which effectively means the dunefoot will get lower for narrow beaches, and the odds of erosion will increase.
            slope0=Dfelev0./Wberm;
            Dfelev=Dfelev0.*min(maxslope./slope0,1);
            slope=Dfelev./Wberm;

            % Offshore wave conditions (or at toe of dynamic profile(?)
            HSo1 = HSi.*sqrt(max(cosd(dPHIi),eps));
            L0 = g/(2*pi)*TPi.^2;
            sqrHsL0=sqrt(HSo1.*L0);
            
            dtdune=dt/dtfractions;
            [~,qs,qss,ql,R,xtill]=dune_erosion(runupform,runupfactor,Wberm,Dfelev,Dcelev,h0,sqrHsL0,SWLi,TPi,Cs,Cstill,xtill,perctill,A,dtdune,dtnum);
            
            %% Aeolian transport rate from beach to dune

            % Critical and shear velocity
            uc=duneAw*sqrt((rhos-rhoa)*g*d50/rhoa);                                     % wind critical shear velocity (m/s)
            uc=uc.*ones(1,n);
            z0=d50/30;                                                                  % Larson 2016
            %z0=0.081*log(DUNE.d50/0.18);                                               % Zingg (1953) Hallin et al. (2019)
            us=k.*uwindi./log(z/z0);                                               % wind shear velocity (m/s)

            % Potential aeolian transport rate (Kg/m/s)
            mov=us>uc;
            mwe(mov)=Kw*sqrt(d50/d50r)*rhoa*us(mov).^2.*(us(mov)-uc(mov))/g;               
            mwe(~mov)=0;

            % Equilibrium transport rate [m3/m/yr] 
            qwe=spyr*mwe/rhos/(1-por);                                                  % adapted by Anouk                                                 

            % Estimation of the dry beach width
            Bdry=max((1-(R+SWLi)./Dfelev),0).*Wberm ;      

            % Fetch length
            cosphiwndloc=cosd(phic-phiwndi);
            FL = max(Bdry./cosphiwndloc,0);

            % Aeolian transport rate
            qw0=qwe.*(1-(exp(-segmaw.*FL))); 
            qw=qw0.*cosphiwndloc; 
            
            % export data
            COAST.qw=(COAST.qw*(tt-1)+qw)/tt;    % compute average qw, i.e. component due to wind transport
            COAST.qs=(COAST.qs*(tt-1)+qs)/tt;    % compute average qs, i.e. dune erosion that is beneficial for the beach 
            COAST.qss=(COAST.qss*(tt-1)+qss)/tt; % compute average qss, i.e. dune erosion that is beneficial for the beach (i.e. the sandy part, discarding the till part)
            COAST.ql=(COAST.ql*(tt-1)+ql)/tt;    % compute average ql, i.e. dune erosion that is not benefecial for the beach 
            COAST.R=max(COAST.R,R);              % use maximum instead of average for the runup
            COAST.SWL=max(COAST.SWL,SWLi);       % use maximum instead of average for the still water level
        end
        
        COAST.qs(TRANSP.shadow)     = 0;      % Set dune transport to zero due to dune-dune shadowing
        COAST.qs(TRANSP.shadow_h)   = 0;      % Set dune transport to zero due to being sheltered behind structures
        COAST.qs(TRANSP.shadowS_hD) = 0;      % Set dune transport to zero due to being sheltered behind revetment
        COAST.qss(TRANSP.shadow)     = 0;      % Set dune transport to zero due to dune-dune shadowing
        COAST.qss(TRANSP.shadow_h)   = 0;      % Set dune transport to zero due to being sheltered behind structures
        COAST.qss(TRANSP.shadowS_hD) = 0;      % Set dune transport to zero due to being sheltered behind revetment
        COAST.ql(TRANSP.shadow)     = 0;      % Set dune transport to zero due to dune-dune shadowing
        COAST.ql(TRANSP.shadow_h)   = 0;      % Set dune transport to zero due to being sheltered behind structures
        COAST.ql(TRANSP.shadowS_hD) = 0;      % Set dune transport to zero due to being sheltered behind revetment   

    end
end

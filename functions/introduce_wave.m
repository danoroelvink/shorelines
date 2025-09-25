function [WAVE]=introduce_wave(WAVE,TIME,COAST,CC)
% function [WAVE]=introduce_wave(WAVE,TIME,COAST,CC)
% 
% Interpolate wave conditions of multiple offshore/nearshore
% wave data sources along the coastline at every timestep.
% If only one wave source is available, or when only static 
% conditions are prescribed then no interpolation is needed.
% 
% INPUT: 
%   WAVE                   : Structure with wave data of the ShorelineS model, of which is used
%        <option 1>
%        .Hs               : Offshore wave height (input = single number)
%        .Tp               : Offshore wave period (input = single number)
%        .PHIw             : Offshore wave direction (input = single number)
%        <option 2>
%        .wvcfile          : filename of wave climate data file (if 'WVC' variable is used)
%        .WVC              : Offshore wave parameters from time-series or wave climate with probabilities
%                            for a wave timeseries with fields .timenum, .Hs, .Tp and .Dir (each with [Nx1] column vector with data)
%                            for a wave climate with fields .Hs, .Tp, .Dir and .Prob (each with [Nx1] column vector with data)
%        .iwc              : index of the wave climate condition that was selected for this timestep 
%                            (it is re-used for interpolation of other elements at the same time instance)
%        <option 3>
%        .c1               : Transport coefficient of the S-Phi curve [m3/yr/°]
%        .c2               : Coefficient for curvature of the S-Phi curve [1/°]
%        .PHIequi          : Equilibrium angle of coastline [°N]
%        .QSoffset         : Tide transport offset value [m3/yr]
%        .PHIf             : Orientation of the lower shoreface [°N]
%   TIME                   : current moment in time 'tnow' and 'tprev' is used (in datnum format [days since 0000-00-01 00:00])
%        .tc               : ratio of adaptive time step that is used (0 means a fixed timestep)
%        .tnow             : current time [in days in datenum format]
%   COAST                  : Structure with x,y locations of coastline points
%   CC                     : Correction values on Hs, PHI accounting for climate change
% 
% OUTPUT:
%   WAVE
%        .HSo              : Significant wave height at considered time instance (m)
%        .PHIo             : Wave direction at considered time instance (radians cartesian)
%        .TP               : Wave period at considered time instance (s)
%        .PHIf0            : Orientation of the lower shoreface (only in case option 2 or 4)
%        .iwc              : index of the wave climate condition that was selected for this timestep 
%                            (which is re-used for interpolation of other elements at the same time instance)
%        .i_mc             : index of currently evaluated coastline element
%        .WVCid0           : grid index of each of the offshore wave stations
%        .spacevaryingwave : index showing whether space varying waves are used (0/1)
%        .c1               : (if RAY-files are used) Offshore wave height (input = single number)
%        .c2               : (if RAY-files are used) Offshore wave period (input = single number)
%        .PHIeq            : (if RAY-files are used) Offshore wave direction (input = single number)
%        .QSoffset         : (if RAY-files are used) Tide transport offset
%        .h0               : (if RAY-files are used) Active height of the profile from the RAY-files [m] 
% 			       
%% Copyright notice	       
%   -------------------------------------------------------------------
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

    %% EXTRACT WAVE CLIMATE DATA
    %    INTERPOLATE WAVES IN TIME
    %    INTERPOLATE WAVES ALONG THE GRID
    %method='alongshore_mapping';   
    method='weighted_distance';   
    if isfield(WAVE,'interpolationmethod') 
        method=WAVE.interpolationmethod; 
    end
    xq=COAST.xq;
    yq=COAST.yq;
    nq=COAST.nq;    
    
    %-------------------------------------------------------------------------------------%
    %% METHOD 1 : SPATIAL AND TEMPORAL INTERPOLATION OF WAVE DATA FILES
    %-------------------------------------------------------------------------------------%
    WVCid0=[]; 
    PHIequi=[]; 
    c1=[]; 
    c2=[]; 
    h0=[]; 
    QSoffset=[]; 
    distw=[];

    % get wave data structure WVC
    WVC=WAVE.WVC;
    
    if ~isempty(WAVE.wvcfile) 
        
        % get wave conditions for particular moment in time (for each of the locations)
        xw=[]; 
        yw=[]; 
        var1.Hs=[]; 
        var1.Tp=[]; 
        var2.Dir=[]; 
 
        if isfield(WVC,'c1')
            var1.c1=[];                   
            var1.c2=[]; 
            var1.h0=[]; 
            var1.QSoffset=[]; 
            var2.PHIequi=[]; 
            %var2.PHIequiQS=[];
        end
        
        fieldnm1=get_fields(var1);
        fieldnm2=get_fields(var2);
        for kk=1:length(WVC) 
            if length(unique(WVC(kk).timenum))>1 
                
                if TIME.tnow>max(WVC(kk).timenum) || TIME.tnow<min(WVC(kk).timenum)
                    dt=WVC(kk).timenum-TIME.tnow;
                    idt=find(mod(dt,365.25)<min(mod(dt,365.25))+1/24,1);
                    dtrewind=dt(idt);
                    if min(mod(dt,365.25))>30 || (max(WVC(kk).timenum)-min(WVC(kk).timenum))<120 % if time-series too short
                        idt=1;
                        dtrewind=dt(1);
                    end
                    if (TIME.tnow+6/24)>min(WVC(kk).timenum) && TIME.tnow<min(WVC(kk).timenum)
                        % extend first and second timepoint to before the first model time instance when it is less than 6 hours
                        WVC(kk).timenum(1)=WVC(kk).timenum(1)-dtrewind;
                        WVC(kk).timenum(2)=WVC(kk).timenum(2)-dtrewind/2;
                        %WAVE.dtrewind(kk)=WAVE.dtrewind(kk)+dtrewind;
                        WAVE.WVC=WVC;                       
                    else
                        % shift the whole time-series by dtrewind in case of larger than 6 hour difference
                        WVC(kk).timenum=WVC(kk).timenum-dtrewind;
                        WAVE.dtrewind(kk)=WAVE.dtrewind(kk)+dtrewind;
                        WAVE.WVC=WVC;
                    end
                end 

                if kk==length(WVC) && sum(WAVE.dtrewind)~=0 && COAST.i_mc==1
                    struct2log('Warning',['   Warning : Current time is not covered by wave time-series. Recycling wave data, using ',datestr(TIME.tnow+WAVE.dtrewind(kk),'yyyy-mm-dd'),' for points ',num2str(find(WAVE.dtrewind),'%1.0f '),'!'],'a');
                end

                % regular wave parameters
                for jj=1:length(fieldnm1)
                    var1.(fieldnm1{jj})(1,kk)=interp1(WVC(kk).timenum,WVC(kk).(fieldnm1{jj}),TIME.tnow); 
                end
              
                % other parameters
                for jj=1:length(fieldnm2)
                    var2.(fieldnm2{jj})(1,kk)=mod(atan2d(interp1(WVC(kk).timenum,sind(WVC(kk).(fieldnm2{jj})),TIME.tnow),...
                                                         interp1(WVC(kk).timenum,cosd(WVC(kk).(fieldnm2{jj})),TIME.tnow)),360);
                end
                WAVE.Prob=1;
            else 
                % wave climate
                if kk==1
                    if COAST.i_mc==1 && TIME.dtsteps==0 
                        if TIME.tc~=0
                            struct2log('Warning','   Warning : Wave climate input should not be combined with a flexible timestep!','a');
                        end
                        if WAVE.mergeconditions==1
                            WAVE.rnd=[];
                            WAVE.iwc=[1:length(WVC(1).Hs)]';
                            WAVE.Prob=WVC(1).Prob(:)./sum(WVC(1).Prob);
                        elseif isfield(WVC(1),'Prob')
                            WAVE.rnd=[];
                            WAVE.iwc=get_randsample(length(WVC(1).Hs),1,WVC(1).Prob);
                            WAVE.Prob=WVC(1).Prob(:)./sum(WVC(1).Prob);
                        else
                            WAVE.rnd=rand;    % random number for drawing from wave climate
                            WAVE.iwc=round((WAVE.rnd*length(WVC(1).Hs)+0.5));
                            WAVE.Prob=ones(length(WVC(1).Hs),1)/length(WVC(1).Hs);
                        end
                    end
                end
                for jj=1:length(fieldnm1)
                    var1.(fieldnm1{jj})(:,kk)=WVC(kk).(fieldnm1{jj})(WAVE.iwc)'; 
                end
                for jj=1:length(fieldnm2)
                    var2.(fieldnm2{jj})(:,kk)=WVC(kk).(fieldnm2{jj})(WAVE.iwc)'; 
                end
            end 
            if isfield(WVC(kk),'x')
               xw(1,kk)=WVC(kk).x;
               yw(1,kk)=WVC(kk).y; 
            end
        end 
        
        % find the right alongshore location for each of the wave climates
        % sort locations alongshore
        % choose only the closest wave climate at each grid cell, to make sure that no multiple climates are enforced on a single grid cell.
        % Make sure that wave climates are not too far from the grid cells. 
        % Otherwise throw out the ones that are at great distance.
        %  - remove wave climate points that are further away (in distance) 
        %  - from the grid point than the distance to another nearby climate point.
        %  - and which are more than 2x the cross-shore distance from the coast than the adjacent climate point
        % interpolate the waves at the right alongshore location    
        if WAVE.spacevaryingwave
            [var1N,var2N,idGRID,distw]=get_interpolation_on_grid(method,xq,yq,xw,yw,var1,var2);
            Hs=var1N.Hs;
            Tp=var1N.Tp;
            PHIw=var2N.Dir;
            WVCid0=idGRID;
            
            if isfield(WVC,'c1')
                c1=var1N.c1';
                c2=var1N.c2';
                h0=var1N.h0';
                QSoffset=var1N.QSoffset';
                PHIequi=var2N.PHIequi';
            end
        else
            Hs=repmat(var1.Hs,[1,length(xq)]); 
            Tp=repmat(var1.Tp,[1,length(xq)]); 
            PHIw=repmat(var2.Dir,[1,length(xq)]); 
            WVCid0=[1:nq]';
        end
        
    %-------------------------------------------------------------------------------------%
    %% METHOD 2 : USE FIXED WAVE HEIGHT AND PULL RANDOM WAVE DIRECTION FROM DIRECTION RANGE
    %-------------------------------------------------------------------------------------%
    else
        if COAST.i_mc==1
        WAVE.rnd=rand;    % random number for drawing from wave climate
        end
        Hs=WAVE.Hso;
        PHIw=WAVE.phiw0+(WAVE.rnd-.5)*WAVE.spread;
        Tp=WAVE.tper;

        WAVE.Hso=WAVE.Hso;        % wave height [m]
        WAVE.phiw0=WAVE.phiw0;    % deep water wave angle in degrees [?N]
        WAVE.spread=WAVE.spread;  % wave spreading [?] (wave_dir from range:  S.phiw0 +/- 0.5*S.spread)
        WAVE.Prob=1;
    end
    
    % store spacevarying wave
    if length(Hs)>1 
        WAVE.spacevaryingwave=1;
    end
    
    if length(PHIw)==1
        PHIw=repmat(PHIw,[1,nq]);
    end
    if length(Hs)==1 
        Hs=repmat(Hs,[1,nq]);
    end
    if length(Tp)==1 
        Tp=repmat(Tp,[1,nq]);
    end
    
    
    %% export variables
    eps=0.01;
    WAVE.HSo=max(Hs,eps)*CC.HScor;
    WAVE.PHIo=mod(PHIw+CC.PHIcor,360.);
    WAVE.TP=max(Tp,eps);
    WAVE.dist=distw;
    WAVE.i_mc=COAST.i_mc;
    try
        WAVE.WVCid0=WVCid0;
    end
    try
        WAVE.PHIeq=mod(PHIequi+CC.PHIcor,360.);   
        WAVE.c1=c1*CC.HScor.^2;
        WAVE.c2=c2;
        WAVE.h0=h0;
        WAVE.QSoffset=QSoffset;
    end
end

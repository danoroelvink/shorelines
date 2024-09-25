function [D]=get_inputfiledata(DATAfile,TIME)
% function [D]=get_inputfiledata(DATAfile,TIME)
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
    
    D=struct;
    D.spacevaryingwave=0;  
    D.timenum=[];
    
    % read content of files
    if ~isempty(DATAfile)
        % re-organize input wave files 
        if ischar(DATAfile)
            DATAfile = {DATAfile};
        end       
        if (size(DATAfile,2)>1 && size(DATAfile,1)==1)
            DATAfile=DATAfile';
        end
        
        %% CONVERT FILES THAT REFERENCE MULTIPLE WAVE CLIMATE FILES/LOCATIONS
        % convert file format if a WVT/WVC file is used as input with reference to a large number of WVT/WVC files 
        if (strcmpi(DATAfile{1}(end-3:end),'.WVT') || strcmpi(DATAfile{1}(end-3:end),'.WVC') ... 
              || strcmpi(DATAfile{1}(end-3:end),'.WVD') || strcmpi(DATAfile{1}(end-3:end),'.WAT') ... 
              || strcmpi(DATAfile{1}(end-3:end),'.WND')) && length(DATAfile)==1
            try 
                Draw=load(DATAfile{1});
            catch 
                fid0=fopen(DATAfile{1},'r');
                DATAfile = {};
                gg=0;
                while ~feof(fid0)
                    line={};
                    try
                        gg=gg+1;
                        line=fgetl(fid0);
                    end
                    info = textscan(line,'%s %f %f');
                    [dirnm,filnm,extnm]=fileparts(info{1}{1});
                    if isempty(dirnm)
                       dirnm='.'; 
                    end    
                    DATAfile{gg,1} = [dirnm,filesep,filnm,extnm];
                    DATAfile{gg,2} = info{2};
                    DATAfile{gg,3} = info{3};
                end
            end
        end
        
        % convert file format if a GKL-file is used as input with reference to a large number of RAY files 
        if strcmpi(DATAfile{1}(end-3:end),'.GKL') && length(DATAfile)==1
            GKLdata=readGKL(DATAfile{1});
            [dirnm,filnm,extnm]=fileparts(DATAfile{1});
            if isempty(dirnm)
                dirnm='.';
            end
            DATAfile = {};
            for gg=1:length(GKLdata.x)
                DATAfile{gg,1} = [dirnm,filesep,GKLdata.ray_file{gg},'.ray'];
                DATAfile{gg,2} = GKLdata.x(gg);
                DATAfile{gg,3} = GKLdata.y(gg);
            end
        end
        
        % read mat-file here 
        readmatfile=1;
        if strcmpi(DATAfile{1}(end-3:end),'.mat') && length(DATAfile)==1
            Draw=load(DATAfile{1});
            if isfield(Draw,'timenum') 
                for gg=2:length(Draw)
                    DATAfile{gg}=DATAfile{1};
                end
                readmatfile=0;
            end
        end
        
        %% LOOP OVER WAVE TIME SERIES POINTS IN SPACE
        kk1=1;
        D(1).spacevaryingwave=size(DATAfile,1)>1;
        for kk=1:size(DATAfile,1)
            fprintf('   reading file : %s\n',DATAfile{kk});
            warning off;
            
            % READ NC-FILES FROM SNAPWAVE
            if strcmpi(DATAfile{kk}(end-2:end),'.nc')
                x=ncread(DATAfile{kk},'station_x');
                y=ncread(DATAfile{kk},'station_y');
                hm0=ncread(DATAfile{kk},'point_hm0');
                tp=ncread(DATAfile{kk},'point_tp');
                wd=ncread(DATAfile{kk},'point_wavdir');
                
                timeunits=ncreadatt(DATAfile{kk},'time','units');
                reftime=datenum(timeunits(15:end)); %reftime=datenum(1979,01,01);
                tt=ncread(DATAfile{kk},'time');
                time=reftime+double(tt)/24/60/60;   % single precision unable to make out seconds->years
                for kk2=1:length(x)
                    D(kk1).x=x(kk2);
                    D(kk1).y=y(kk2);
                    D(kk1).timenum=time(:);
                    if length(x)==size(hm0,1)
                        D(kk1).Hs=squeeze(hm0(kk2,:))';
                        D(kk1).Tp=squeeze(tp(kk2,:))';
                        D(kk1).Dir=squeeze(wd(kk2,:))';
                    else
                        D(kk1).Hs=squeeze(hm0(:,kk2));
                        D(kk1).Tp=squeeze(tp(:,kk2));
                        D(kk1).Dir=squeeze(wd(:,kk2));
                    end
                    kk1=kk1+1;
                end
                D(1).spacevaryingwave=length(D)>1;
            
            % READ WAVE TIME-SERIES OR CLIMATE DATA (FOR RUNUP)
            elseif strcmpi(DATAfile{kk}(end-3:end),'.WVC') || strcmpi(DATAfile{kk}(end-3:end),'.WVT') ...
                   || strcmpi(DATAfile{kk}(end-3:end),'.WVD') || strcmpi(DATAfile{kk}(end-3:end),'.TXT')
                Draw=load(DATAfile{kk});
                if strcmpi(DATAfile{kk}(end-3:end),'.WVT') || Draw(1,1)>1000 
                    % offshore wave time-series
                    if length(num2str(Draw(1,1)))==12
                        D(kk).timenum=datenum([num2str(Draw(:,1))],'yyyymmddHHMM');
                    elseif length(num2str(Draw(1,1)))<12
                        D(kk).timenum=datenum([num2str(Draw(:,1))],'yyyymmdd');
                    else
                        D(kk).timenum=datenum([num2str(Draw(:,1))],'yyyymmddHHMMSS');       
                    end
                    D(kk).Hs=interpNANs(Draw(:,2));
                    D(kk).Tp=interpNANs(Draw(:,3));
                    D(kk).Dir=interpNANsDIR(Draw(:,4));
                else % strcmpi(DATAfile{kk}(end-3:end),'.WVC') 
                    % offshore wave climate
                    D(kk).timenum=[];
                    D(kk).Hs=interpNANs(Draw(:,1))';
                    D(kk).Tp=interpNANs(Draw(:,2))';
                    D(kk).Dir=interpNANsDIR(Draw(:,3))';
                    if size(Draw,2)>3
                        D(kk).Prob=interpNANs(Draw(:,4))';
                        
                        % scale to 1.0
                        if sum(D(kk).Prob)>1.1 && sum(D(kk).Prob)<100      % percentage -> fraction of 1
                            D(kk).Prob=D(kk).Prob/100.; 
                        elseif sum(D(kk).Prob)>200                      % days -> fraction of 1
                            D(kk).Prob=D(kk).Prob/365.; 
                        else                                          % other -> fraction of 1
                            D(kk).Prob=D(kk).Prob/sum(D(kk).Prob);
                        end
                    end
                end
            
            % READ WATER-LEVEL DATA
            elseif strcmpi(DATAfile{kk}(end-3:end),'.WAT') 
                Draw=load(DATAfile{kk});
                if Draw(1,1)>1000
                    % read waterlevel time-series data
                    if length(num2str(Draw(1,1)))==12
                        D(kk).timenum=datenum([num2str(Draw(:,1))],'yyyymmddHHMM');
                    elseif length(num2str(Draw(1,1)))<12
                        D(kk).timenum=datenum([num2str(Draw(:,1))],'yyyymmdd');
                    else
                        D(kk).timenum=datenum([num2str(Draw(:,1))],'yyyymmddHHMMSS');       
                    end
                    D(kk).swl=interpNANs(Draw(:,2));
                else
                    % read waterlevel climate data
                    D(kk).timenum=[];
                    D(kk).swl=interpNANs(Draw(:,1))';
                    if size(Draw,2)>2
                        D(kk).Prob=interpNANs(Draw(:,2))';    
                        
                        % scale to 1.0
                        if sum(D(kk).Prob)>1.1 && sum(D(kk).Prob)<100      % percentage -> fraction of 1
                            D(kk).Prob=D(kk).Prob/100.; 
                        elseif sum(D(kk).Prob)>200                      % days -> fraction of 1
                            D(kk).Prob=D(kk).Prob/365.; 
                        else                                          % other -> fraction of 1
                            D(kk).Prob=D(kk).Prob/sum(D(kk).Prob);
                        end
                    end
                end
            
            % READ WIND DATA
            elseif strcmpi(DATAfile{kk}(end-3:end),'.WND') || strcmpi(DATAfile{kk}(end-3:end),'.WDT') || strcmpi(DATAfile{kk}(end-3:end),'.WDC')  
                Draw=load(DATAfile{kk});
                if Draw(1,1)>1000
                    % wind time-series data
                    if length(num2str(Draw(1,1)))==12
                        D(kk).timenum=datenum([num2str(Draw(:,1))],'yyyymmddHHMM');
                    elseif length(num2str(Draw(1,1)))<12
                        D(kk).timenum=datenum([num2str(Draw(:,1))],'yyyymmdd');
                    else
                        D(kk).timenum=datenum([num2str(Draw(:,1))],'yyyymmddHHMMSS');       
                    end
                    D(kk).uz=interpNANs(Draw(:,2));
                    D(kk).Dir=interpNANs(Draw(:,3));
                else
                    % wind climate data
                    D(kk).timenum=[];
                    D(kk).uz=interpNANs(Draw(:,1))';
                    D(kk).Dir=interpNANsDIR(Draw(:,2))';
                    if size(Draw,2)>2
                        D(kk).Prob=interpNANs(Draw(:,3))';    
                        
                        % scale to 1.0
                        if sum(D(kk).Prob)>1.1 && sum(D(kk).Prob)<100      % percentage -> fraction of 1
                            D(kk).Prob=D(kk).Prob/100.; 
                        elseif sum(D(kk).Prob)>200                      % days -> fraction of 1
                            D(kk).Prob=D(kk).Prob/365.; 
                        else                                          % other -> fraction of 1
                            D(kk).Prob=D(kk).Prob/sum(D(kk).Prob);
                        end
                    end
                end
            
            % TIME-SERIES OR STATIC RAY FILES
            elseif strcmpi(DATAfile{kk}(end-3:end),'.RAY')
                RAYdata=readRAY(DATAfile{kk,1});
                D(kk).time     = [];  
                D(kk).timenum  = TIME.timenum0;    % use as fall-back option the starttime of the simulation as starttime of the datafile
                try  
                    D(kk).time = RAYdata.time;     % import/read time (only in case of time-series ray-file)
                end
                if ~isempty(RAYdata.time)
                    D(kk).timenum  = D(kk).timenum+RAYdata.time/24;
                end
                D(kk).PHIequi  = RAYdata.Cequi;                  % angle of equilibrium also accounting for the net currents which are in QSoffset (=computeEQUI(D(kk).equi,D(kk).c1,D(kk).c2,D(kk).hoek,D(kk).QSoffset) 
                D(kk).c1       = RAYdata.c1*10^6;     
                D(kk).c2       = RAYdata.c2; 
                D(kk).h0       = RAYdata.h0;    
                D(kk).fshape   = RAYdata.fshape; 
                D(kk).PHIf     = RAYdata.hoek;
                D(kk).QSoffset = zeros(size(D(kk).PHIequi));
                if ~isempty(RAYdata.QSoffset)
                    D(kk).QSoffset = RAYdata.QSoffset; % timeseries ray with tide offset
                end
                D(kk).Hs        = (abs(D(kk).c1).^0.5)/200.0;    % use a proxy     Hs=sqrt(c1)/200
                D(kk).Tp        = D(kk).Hs*2+3;                  % use a proxy     Tp=Hs*2+3
                D(kk).Dir       = RAYdata.hoek-RAYdata.equi;       % angle of equilibrium accounting only for the wave part (i.e. hoek-equi)              % use a proxy     Dir=Cequi
            
            % TIME-SERIES FROM WAVE-MAT-FILE FOR DATA ASSIMILATION OR COUPLING
            elseif strcmpi(DATAfile{kk}(end-3:end),'.mat')
                if readmatfile==1
                    % read Draw for each location, if not read earlier
                    Draw=load(DATAfile{kk});
                end
                if isfield(Draw,'WaveDA_record')
                    D(kk).timenum=Draw.WaveDA_record(:,1);
                    D(kk).Hs=interpNANs(Draw.WaveDA_record(:,2));
                    D(kk).Tp=interpNANs(Draw.WaveDA_record(:,3));
                    D(kk).Dir=interpNANsDIR(Draw.WaveDA_record(:,4));
                elseif isfield(Draw,'WVCn')
                    D(kk).timenum=datenum([num2str(Draw.WVCn(:,1))],'yyyymmdd');
                    D(kk).Hs=interpNANs(Draw.WVCn(:,2));
                    D(kk).Tp=interpNANs(Draw.WVCn(:,3));
                    D(kk).Dir=interpNANsDIR(Draw.WVCn(:,4));
                elseif isfield(Draw,'timenum')
                    % A mat-file with a matlab-structure with wavedata-fields for multiple locations
                    D(kk).timenum=interpNANs(Draw(kk).timenum(:));
                    D(kk).Hs=interpNANs(Draw(kk).Hs(:));
                    D(kk).Tp=interpNANs(Draw(kk).Tp(:));
                    D(kk).Dir=interpNANsDIR(Draw(kk).Dir(:));
                    D(kk).x=interpNANsDIR(Draw(kk).x);
                    D(kk).y=interpNANsDIR(Draw(kk).y);
                end
            
            % IF NO SUITBALE INPUT IS PROVIDED
            else
                fprintf(['   Warning : Wave condition file ''',DATAfile{kk},''' cannot be read. Please check the formatting! \n'])
            end
            
            % ADD LOCATION X AND Y
            if ~strcmpi(DATAfile{kk}(end-2:end),'.nc')
                D(kk).x = kk;
                D(kk).y = kk;
                if size(DATAfile,2)==3
                    D(kk).x = DATAfile{kk,2};
                    D(kk).y = DATAfile{kk,3};
                else
                    if kk>1
                        fprintf([' Warning : please specify an X and Y location for ''',DATAfile{kk},''' (i.e. second and third column of S.wvcfile) \n']);
                    end
                end
            end
        end 
        
        %% check if timeseries covers valid time range (after model start / 'timenum0')
        for kk=1:length(D)
            if ~isempty(D(kk).timenum)
                if D(kk).timenum(1) > TIME.timenum0
                    fprintf([' Warning : The model start time is not covered by the wave data for timeseries file ''',DATAfile{kk},'''\n']);
                end
            end
        end
        
        %% determine additional counter which tracks the index of the next timestep after the current time. So, the index at t0 +1.
        % this is used for determining the adaptive time step
        % indxw should always be the index of the time+dt (i.e. timeseries timestep that comes after the current time).
        for kk=1:length(D)
            D(kk).indxw=find(D(kk).timenum>TIME.timenum0,1); 
            D(kk).indxw=min(D(kk).indxw,length(D(kk).timenum));
        end
    end
end

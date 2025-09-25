function [S]=initialize_defaultvalues(S0)
% function [S]=initialize_defaultvalues(S0)
%
% INPUT:
%    S       Default settings of the ShorelineS model (as specified in this function)
%    S0      User input given to the ShorelineS model
%
% OUTPUT:
%    S       Combined settings with defaults replaced by user settings
%
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2020 IHE Delft & Deltares
%
%       Dano Roelvink, Ahmed Elghandour
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

    fprintf('%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');  
    fprintf('%s\n','%                  RUN SHORELINES                   %');  
    fprintf('%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');  
    fprintf('%s\n','  Initialize default values');

    if (isoctave)
        warning('off','all');
    end
    
    %% ---------------------- simulation wave parameters ----------------------        
    S.hso=1;                                                                   % wave height [m]
    S.tper=6;                                                                  % peak wave period [s]
    S.phiw0=330;                                                               % deep water wave angle in degrees [?N]
    S.spread=90;                                                               % wave spreading of the wave climate [�] (i.e. variation over time), affecting the wave_dir from range S.phiw0 +/- 0.5*S.spread
    S.dirspr=12;                                                               % directional spreading of the wave [�] (i.e. of the 2D wave spectrum at a moment in time) affecting the diffraction at structures
    S.wvcfile='';                                                              % wave time-series <-leave empty to use wave parameters ('S.hso', 'S.phiw0' and 'S.spread')
    S.ddeep=25;                                                                % Waterdepth the location of wave climate, corresponding with S.hso [m]
    S.dnearshore=8;                                                            % Waterdepth at the 'dynamic boundary', corresponding with S.phif [m]
    S.randomseed=-1;                                                           % Seed to generate randomized series
    S.randomseedsettings='';                                                   % Random generator settings
    S.interpolationmethod='weighted_distance';  %'alongshore_mapping';         % Method for interpolating wave data of multiple wave stations on the coast (either 'alongshore_mapping' / 'weighted_distance')
    S.mergeconditions=0;                                                       % switch to turn on the use of all wave, wind and runup climate conditions simultaneously (of a WVC-file/WND-file/WVD-file/WLC-file) instead of just using a single wvc-condition per timestep (no effect for wvt-files)
    %% ---------------------simulation tide parameters ------------------------
    S.tidefile='';                                                             % tidal forcing needed as input for tidemodule S.trform = 'tideprof'. Input colums [xstat ystat etaM2 etaM4 detadsM2 detadsM4 phiM2 phiM4 kM2 kM4 ss]
    S.tideprofile='';                                                          % cross-shore profile needed as input for tidemodule S.trform = 'tideprof'
    S.tidedx=10;                                                               % cross-shore distance resolution needed as input for tidemodule S.trform = 'tideprof'
    S.tiden=[];                                                                % manning coefficient for wave induced currents in the tide module [s/m^(1/3)] tidemodule S.trform = 'tideprof'
    %% ------------------- simulation coastline definition --------------------        
    S.ldbcoastline='';                                                         % file with initial coastline shape ([Nx2] ASCII FILE WITHOUT HEADER), or directly [Nx2] data <- leave empty to use interactive mode!
    S.xmc='';                                                                  % x-coordinates of coastline [Nx1], with NaNs separating the elements <- leave empty to use interactive mode!
    S.ymc='';                                                                  % y-coordinates of coastline [Nx1], with NaNs separating the elements <- leave empty to use interactive mode!
    S.gisconvention=0;                                                         % switch to swap coastline definition according to GIS convention (i.e. counter clockwise)
    S.ds0=100;                                                                 % initial space step [m]
    S.griddingmethod=2;                                                        % method for regenerating the grid (1: only splitting and merging cells if they are too small, 2: uniform grid regeneration if criteria for gridsize or exceeded)
    S.d=10;                                                                    % active profile height [m]
    S.phif=[];                                                                 % Orientation of the foreshore [?N] (in degrees) <- only relevant for 'KAMP', 'MILH' or 'VR14'
    S.maxangle=60;                                                             % maximum coastline re-orientation between individual grid cells (affecting spit width and stabilizing small scale features in case of dense grids)
    S.preserveorientation=0;                                                   % switch for preserving the original orientation of the shoreline (either 0 or 1, by default this is switched off)
    %% ------------------- simulation transport parameters --------------------    
    S.trform='CERC';                                                           % switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
    S.b=1e6;                                                                   % CERC : coeff in simple cerc formula
    S.qscal=1;                                                                 % calibration factor of the transport (works for all transport formulas)
    S.d50=2.0e-4;                                                              % median grain diameter [m] (KAMP & MILH & VR14)
    S.d90=3.0e-4;                                                              % 90th percentile grain diameter [m] (KAMP & MILH & VR14)
    S.porosity=0.4;                                                            % S.porosity (typically 0.4) [-] (KAMP & MILH & VR14)
    S.tanbeta=0.03;                                                            % mean bed slope [ratio 1/slope] (KAMP & MILH & VR14)
    S.tanbetasetup=1;                                                          % bed slope value, scaling the effect of the water-level setup driven currents impact on sediment transport (dHs/ds), with 1 as default for a very small impact on transport (for a steep slope)
    S.rhos=2650;                                                               % density of sand [kg/m3] (KAMP & MILH & VR14)
    S.rhow=1025;                                                               % density of water [kg/m3] (KAMP & MILH & VR14)
    S.g=9.81;                                                                  % gravitational acceleration [m2/s] (KAMP & MILH & VR14)
    S.alpha=1.8;                                                               % calibration factor for point of braking (S.alpha = 1.8 for Egmond data) (KAMP & MILH & VR14)
    S.gamma=0.72;                                                              % breaking coefficient (Hs/h) with 5% breaking waves (KAMP & MILH & VR14)
    S.pswell=20;                                                               % Percentage swell (between 0 - 100) [-] (VR14)
    S.ks=0.05;                                                                 % roughness height [m] 
    S.hclosure=8;                                                              % depth-of-closure [m] 
    S.cf=0.0023;                                                               % roughness factor [-] used in transport formulation TIDEPROF
    S.n=0.02;                                                                  % manning roughness coefficient [s/m^(1/3)], if specified it is used instead of the 'cf' in the tide module, then cf=9.81*(n^2)./h.^(1/3)
    S.acal=0.2;                                                                % calibration coefficient for Soulsby Van Rijn [-]
    S.hmin=0.1;                                                                % minimum depth [m] used in transport formulation TIDEPROF
    S.sphimax=[];                                                              % the computation of the maximum angle of the waves can be bypassed by specifying this sphimax value (e.g. at sphimax=42) or by using 'auto' which means that it is computed only at t0
    S.relaxationlength=[];                                                     % length over which transport decelerates in meters, which adds inertia to the longshore current. It scales linearly with the wave height below 1m waves.
    S.suppresshighangle=0;                                                     % switch 0/1 to disable the high-angle instabilities by limiting the transport angle to the critical high-angle orientation (when it is set at 1)
    %% ------------------- simulation time steps & numerical-------------------
    S.tc=1;                                                                    % switch for using adaptive time step (0/1). The tc-value gives a fraction of the automatically computed timestep. So when S.tc=0.5 is used, it will use a timestep which is half of the maximum step that is permissable. S.tc=0 means that the fixed step S.dt is used instead.
    S.dt=0;                                                                    % value of a fixed timestep for the coastline processes [fraction of a year]. It is only used when S.tc=0;
    S.dtdune=[];                                                               % value of a fixed timestep for dune process [fraction of a year]. It should always be less than 'dt' of the coastline;
    S.reftime='2020-01-01';                                                    % reference/start time of the simulation ('yyyy-mm-dd') <- leave empty to use t=0
    S.endofsimulation='2040-01-01';                                            % the end time of the simulation ('yyyy-mm-dd')
    S.twopoints=1;                                                             % switch for 'S.twopoints approach' which determines the type of response during high-angle wave events (default=1)
    S.smoothfac=0;                                                             % smoothing factor used to re-arrange grid every timestep (only for griddingmethod==1), with reasonable values between 0 and 0.1.
    S.smoothrefrac=0;                                                          % smoothing fraction of the coastline orientation (PHIcs) which is used only for the refraction of waves from nearshore location (TDP) to point of breaking (BR), and not for the transport computation. Values between 0 and 1.
    %% -------------------------- boundary condition -------------------------- 
    S.boundaryconditionstart='Fixed';                                          % left boundary condition (either 'Fixed', 'Closed', 'Angleconstant' or 'Periodic') for non-cyclical elements. An additional argument with the specified coastangle or transport rate can be given by providing a cell-format. For example, {'Closed',-6000} to specify a negative transport of -6000 m3/yr or {'Angleconstant',300} to enforce a coastangle of 300캮. 
    S.boundaryconditionend='Fixed';                                            % right boundary condition (either 'Fixed', 'Closed', 'Angleconstant' or 'Periodic') for non-cyclical elements. An additional argument with the specified coastangle or transport rate can be given by providing a cell-format. For example, {'Closed',-6000} to specify a negative transport of -6000 m3/yr or {'Angleconstant',300} to enforce a coastangle of 300캮. 
    S.cyclic=0;                                                                % switch for cyclic boundary conditions (e.g. for mud)
    %% ------------------------ climate change impact -------------------------
    S.ccslr=[];                                                                % climate impacted rise in sea level (SLR). This can be a constant rate of sea level rise (e.g. S.ccslr=0.002 m/yr) or a table with the absolute sea level against time [Nx2]. Wherein 'time in datenum format' and 'sea level with respect to initial situation'). Note that the rates per year are computed automatically! For example, the rate of SLR goes up from 2mm/yr to 10mm/yr in this table: S.ccslr=[t1, 0.002; t2,0.003; t3,0.005; t4,0.01]; The S.tanbeta is used as 'slope angle' for the BRUUN rule. 
    S.cchs=[];                                                                 % climate impacted increase of the wave height (HS). This can be a constant rate per year (e.g. S.cchs=0.001, which is +0.1% increase in HS per year) or a table with for various time instances the relative change in wave height with respect to the initial situation [Nx2]. Wherein 'time in datenum format' and 'relative change in wave height w.r.t. initial situation' as a factor. For example, the HS may increase by up to 40% up till 't4' in this table S.cchs=[t1, 0.0; t2,0.1; t3,0.25; t4,0.4]; 
    S.ccdir=[];                                                                % climate impacted change in wave direction (DIR). This can be a constant rate (e.g. S.ccdir=0.05 째/yr) or a table with the wave direction change with respect to the initial situation against time [Nx2]. Wherein 'time in datenum format' and 'relative change in wave direction w.r.t. initial situation' as a # degrees. For example, the wave direction may be adjsuted by up to 0째 up till +3째 at 't4' in this table S.ccdir=[t1, 0째; t2,+0.5째; t3,+1.5째; t4,+3째]; 
    %% ----------------------------- structures -------------------------------
    S.struct=1;                                                                % switch for using hard structures (0/1)
    S.ldbstructures='';                                                        % file with hard structures (ASCII FILE WITHOUT HEADER with [Nx2]) / or [Nx2] matrix <- leave empty to use interactive mode!
    S.xhard=[];                                                                % x-coordinates of hard structures ([1xN]) seperated by nans <- leave empty to use interactive mode!
    S.yhard=[];                                                                % y-coordinates of hard structures ([1xN]) seperated by nans <- leave empty to use interactive mode!
    S.structtype={};                                                           % type of structures
    S.aw=1.27;                                                                 % factor for determining depth of closure at bypassing groyne (1.27 if time series is used) This value is used by default.
    S.awfixedhs=5;                                                             % factor for determining depth of closure at bypassing groyne for a representative Hs is used instead of a climate or timeseries. This value is used instead of 'aw' if S.wvcfile is empty.
    S.bypasscontrfac=1;                                                        % scaling factor for the transport bypass at a groyne as a result of contraction of the flow (always >=1). Setting the bypass contraction factor larger than 1 means that the accretion does not go to the tip of the structure. 
    S.bypassdistpwr=1;                                                         % this power controls the distribution of bypassed sediment of groynes in the shadow zone. It will be a triangle with pwr=1 (default), and more sediment will end up closer to the structure for a higher power (e.g. pwr=2).
    S.submerged=0;                                                             % switch (0/1) to use submerged groynes (experimental feature)
    S.groinelev=[];                                                            % elevation of submerged groynes
    %% ------------- wave transmission over submerged breakwater --------------
    S.transmission=0;                                                          % switch to allow for wave transmission over breakwaters (0/1)
    S.transmform='angr';                                                       % wave transmission approaches for rough structures: d'Angremond (1996) 'angr', Van Gent et al. for rough structures (2023) 'gent' or Seabrook and Hall (1998) 'seabrhall'        
    S.transmdir=1;                                                             % switch to adapt wave direction after transmission over breakwater scaling with percentage transmission [S.transmdir=1], or keep original wave direction [S.transmdir=0];
    S.transmfile='';                                                           % file with characteristics of submergable breakwaters ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
    S.transmbwdepth=[];                                                        % depth at which breakwater is positioned [m] Positive downward.
    S.transmcrestheight=[];                                                    % breakwater height [m] w.r.t. MSL Positive upward, above MSL.
    S.transmslope=[];                                                          % breakwater slope [-] 
    S.transmcrestwidth=[];                                                     % breakwater crest width [m] 
    S.transmd50=1;                                                             % D50 of armour layer, needed as input for Seabrook and Hall formulation 'SeabrHall'. Put to 1 m as default.                
    %% ------------------------- permeable structures -------------------------
    S.perm=0;                                                                  % switch for using permeable structures (0/1)
    S.ldbpermeable='';                                                         % file with permeable structures ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
    S.xperm=[];                                                                % x-coordinates of permeable structures ([1xN]) seperated by nans <- leave empty to use interactive mode!
    S.yperm=[];                                                                % y-coordinates of permeable structures ([1xN]) seperated by nans <- leave empty to use interactive mode!
    S.wavetransm=[1];                                                          % wave transmission at permeable structures
    %% ----------------------------- revetments -------------------------------
    S.revet=1;                                                                 % switch for using revetments (0/1)
    S.ldbrevetments='';                                                        % LDB with hard structures ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
    S.xrevet=[];                                                               % x-coordinates of revetments ([1xN]) seperated by nans <- leave empty to use interactive mode!
    S.yrevet=[];                                                               % y-coordinates of revetments ([1xN]) seperated by nans <- leave empty to use interactive mode!
    S.iterrev=5;                                                               % number of iterations used to determine bypassing along a revetment (sometimes with very long revetments, with spatially varying conditions a large value is needed for a stable computation)
    S.critwidth=5;                                                             % critical width used for determing coastal width in front of revetments. This scales the bypass accordingly!
    %% --------------------------- wave diffraction ---------------------------
    S.diffraction=0;                                                           % switch for using wave diffraction (0/1)
    S.wdform='Roelvink';                                                       % wave diffraction approach to treat angles (Roelvink(default), Hurst). With Hurst the omegat factor is set to 0 and the rotfac to 1, and directional spreading is handled with discrete directional sectors later in the wave_difffraction_coeff function. 
    S.kdform='Kamphuis';                                                       % computation of kd according to 'Kamphuis' or 'Roelvink' analytical approx. The Roelvink method has been derived on the basis of our data. The Kamphuis method gives a smoother result for the wave energy distribution. 
    S.rotfac=1.5;                                                              % factor determining the rotation of waves due to diffraction
    S.diffdist=[];                                                             % specify the maximum distance of structures to coastline points with diffraction (limits influence area)
    S.diffsmooth=0;                                                            % smooth the computed diffraction (re-orientation) [0 - 1] (0 means no smoothing, 1 means smoothed once)
    %% ---------------------------- nourishments ------------------------------
    S.nourish=0;                                                               % switch for using nourishments (0, 1 or 2) -> 2 is common method used
    S.growth=1;                                                                % calibration of nourishment growth rate
    S.norfile='';                                                              % file with nourishments (recommended method), which stores a table with at each line the properties of a nourishment (x1,y1,x2,y2,t1,t2,volume)
    S.nourmethod='default';                                                    % The 'default' method uses only the begin and end point of the nourishment to identify where the nourishment needs to take place, while the 'complex' method sub-diveds the nourishment in 20 smaller parts which eahc can be attributed to a specific coastal element.
    S.ldbnourish='';                                                           % LDB with nourishment locations ([Nx2] ASCII FILE WITHOUT HEADER) (i.e. polygon around relevant grid cells) <- leave empty to use interactive mode!
    S.nourratefile='';                                                         % nourishment rates placed in order for each nourishment polygons (not used if .NOR file is available
    S.nourstartfile='';                                                        % nourishment start dates placed in order for each nourishment polygons (not used if .NOR file is available
    S.nourendfile='';                                                          % nourishment end dates placed in order for each nourishment polygons (not used if .NOR file is available
    S.nourrate=100;                                                            % rate of nourishing (not used if .NOR file is available)
    %% ------------------------ shoreface nourishments ------------------------
    S.fnourish=0;                                                              % switch for shoreface nourishments (0/1)
    S.fnorfile='';                                                             % shoreface nourishment input file 
    S.sal=35;                                                                  % salinity of the water [ppm salt]
    S.temp=10;                                                                 % temperature of the water [degrees celsius]
    S.mb=-0.5;                                                                 % coefficient 1 for shoreface nourishments supply over time
    S.labda0=0.56E-6;                                                          % coefficient 2 for shoreface nourishments supply over time
    S.k=[];                                                                    % coefficient 3 for shoreface nourishments supply over time
    %% -------------------------- sources and Sinks ---------------------------
    S.sourcessinks='';                                                         % sources and sinks rates
    S.ssfile='';                                                               % sources and sinks file
    %% -------------------- aeolian transport to the dunes --------------------
    S.dune=0;                                                                  % switch for using estimate dune evolution (0/1)
    S.ldbdune='';                                                              % filename with dune parameters (.dun file) [Nx5] with xdune,ydune,wberm,dfelev,dcelev in the columns (ASCII FILE WITHOUT HEADER)  <- leave empty to use default values 
    S.kf=0.02;                                                                 % friction coefficient
    S.cs=5e-4;                                                                 % erosion coefficient of the dunes during storms, which scales the rate of erosion
    S.cstill=5e-6;                                                             % erosion coefficient of the dunes during storms for dunes with consolidated till layers, which scales the rate of erosion
    S.xtill=[];                                                                % thickness of layer with sand at the seaward side of the dunes (which can be eroded first) before a till layer is exposed. Leave empty (default) if dunes are 100% sandy (i.e. no till layers), or use a very large thickness then. A spatially varying value can be specified by using a [Nx3] table with x,y,xtill per location at each line. 
    S.perctill=80;                                                             % percentage of till that is present in the till layer of the dunes (i.e. when dune erosion is beyond xtill). A spatially varying value can be specified by using a [Nx3] table with x,y,perctill per location at each line. 
    S.d50r=2.5e-4;                                                             % median reference grain size
    S.rhoa=1.225;                                                              % air density []
    S.duneaw=0.1;                                                              % coefficient (Bagnold, 1937). A larger value increases the threshhold for aeolian transport, and therefore decreases supply
    S.kw=4.2;                                                                  % empirical coefficient (Sherman et al. 2013)
    S.k=0.41;                                                                  % von Karman's coefficient
    S.segmaw=0.1;                                                              % empirical factor used for scaling impact of the fetch length
    S.maxslope=1/15;                                                           % the maximum slope angle (1:slope) with default of 1:20. The dunefoot height is lowered if the beach gets too steep (preserving the max slope).
    S.aoverwash=3;                                                             % rate of overwash
    S.xdune=0;                                                                 % x-coordinates where wberm and dfelev are given (only used if ldbdune is empty)
    S.ydune=0;                                                                 % y-coordinates where wberm and dfelev are given (only used if ldbdune is empty)
    S.wberm=50;                                                                % initial condition for the berm width (m), relevant for the fetch of the eolian transport & the run-up of Stockdon and Ghonim (only used if ldbdune is empty)
    S.dfelev=3;                                                                % dune foot elevation w.r.t. MSL (m) (only used if ldbdune is empty)
    S.dcelev=8;                                                                % dune crest elevation w.r.t. MSL (m) (only used if ldbdune is empty)
    S.csmodel='';                                                              % file for the input of the CS-model (if non-empty, then the CS-model is used for the dunes isntead of the regular dune model)
    %% --------------------------- wind conditions ----------------------------
    S.wndfile='';                                                              % wind time-series filename <-leave empty to use 'S.uz', 'S.phiwnd0' and 'S.spread'
    S.cd=0.002;                                                                % wind drag coefficient (-)
    S.uz=8;                                                                    % wind velocity at z (m)
    S.z=10;                                                                    % elevation of measured wind data
    S.phiwnd0=330;                                                             % wind angle [degN]
    %% --------------------- still water levels and run-up---------------------
    S.runupform='Stockdon';                                                    % switch for run-up formulation (default 'Stockdon', or 'Ghonim', or 'Larson')
    S.runupfactor=1;                                                           % tuning factor for runup (linear)
    S.watfile='';                                                              % Water levels time series file with [Nx2] date/time in 'yyyymmddHHMM' and waterlevel relative to MSL
    S.wvdfile='';                                                              % Wave height time series file with [Nx2] date/time in 'yyyymmddHHMM' and wave height, wave period and wave direction
    S.swl0=0;                                                                  % Fixed still water level relative to MSL. This value will be added to any water-level data provided in a WAT-file (i.e. values in a .WLC or WLT file)
    %%------------------------- Sediment limitations --------------------------
    S.sedlim=1;                                                                % switch for sediment limiter (0/1)
    S.ldbsedlim='';                                                            % file with sediment limitation coordinates and properties, with [Nx3] specification. The columns have the x and y coordinates and width at which transport starts to reduce gradually to zero, separated with NaN's for different sections.
    S.xsedlim=[];                                                              % alternative input with specifying x-coordinates of areas with sediment limitation separately, as [Nx1] separated with NaN's for different sections
    S.ysedlim=[];                                                              % alternative input with specifying y-coordinates of areas with sediment limitation separately, as [Nx1] separated with NaN's for different sections
    S.widthsedlim=[];                                                          % alternative input with specifying width at which transport starts to reduce for areas with sediment limitation separately, as [Nx1] separated with NaN's for different sections
    %% --------------------------- mud properties -----------------------------
    S.mud=0;                                                                   % mud transport option (0/1)
    S.mudtaucr=0.3;                                                            % critical shear stress for erosion (N/m2)
    S.mudm=1.e-4;                                                              % erosion rate (kg/m2/s)
    S.mudb=1000;                                                               % muddy transport zone width (m)
    S.mudmhw=1;                                                                % mhw level [m]
    S.mudmsl=0;                                                                % msl level [m]
    S.mudtfm=10;                                                               % tfm value
    S.mudw=1e-4;                                                               % fall velocity of the muddy sediment [m/s] 
    S.mudbfcrit=500;                                                           % critical mudflat width [m]
    S.mudbmmin=1;                                                              % minimum for the width of the muddy transport zone [m]
    S.mudbmmax=100000;                                                         % maximum for the width of the muddy transport zone [m]
    S.ldbriverdisch='';                                                        % river discharge file with .riv extension, containing [xriv1, yriv1, xriv2, yriv2, tstart, tend, rate]
    S.ldbmangrove='';                                                          % mangrove definition file with .mgv extension, containing [xmgv, ymgv, Bf, Bm, Bfm] 
    %% ----------------------- physics of spit width --------------------------
    S.spitmethod='default';                                                    % Overwash method
    S.spitwidth=50;                                                            % width of tip of spit (used for overwash)
    S.spitheadwidth=200;                                                       % width of tip of spit (used for upwind correction)
    S.owscale=0.1;                                                             % scales the rate of the overwash per timestep (i.e. what part of the deficit is moved to the backbarrier)
    S.owtimescale=0.0;                                                         % timescale for overwash (i.e. what part of the deficit is moved to the backbarrier)
    S.spitdsf=S.d*0.8;                                                         % underwater part of active height for shoreface -> used only in spit-width function
    S.spitdbb=0.5*S.spitdsf;                                                   % underwater part of active height for back-barrier -> used only in spit-width function
    S.bheight=2;                                                               % berm height used for overwash funciton (i.e. added to Dsf or Dbb)
    S.tideinteraction=false;                                                   % tide interaction switch (0/1)
    S.waveinteraction=false;                                                   % wave interaction switch (0/1)
    S.wavefile='';                                                             % wave table (.mat)
    S.surfwidth=1000;                                                          % width of surf zone, where to update the bathymetry
    S.bathyupdate='';                                                          % the dates when the bathymetry should be updated, the input should be in dates form, can accept more than one  {'yyyy-mm-dd'};
    %% ------------------------------- channel --------------------------------
    S.channel=0;                                                               % switch for migrating inlet (0/1)
    S.channelwidth=550;                                                        % target channel width
    S.channelfac=0.08;                                                         % adaptation factor, scales the response of the channel
    S.channeldischrate=0;                                                      % discharge rate
    S.channeldischr=300;                                                       % discharge rate of river
    S.ldbchannel=[];                                                           % file with initial channel axis x and y coordinates, or directly a [Nx2] matrix (option 1)
    S.rero=300;                                                                % radius of flood delta influence on coast
    S.rdepo=600;                                                               % radius of inlet influence on flood delta
    S.tscale=1;                                                                % timescale of flood delta filling
    S.xrmc='';                                                                 % x-coordinates of rivers (option 2)
    S.yrmc='';                                                                 % y-coordinates of rivers (option 2)
    %% ----------------------------- flood delta ------------------------------
    S.flooddelta=0;                                                            % switch (0/1) for flood delta
    S.ldbflood=[];                                                             % wide outline of potential flood delta deposits
    S.xfloodpol=[];                                                            % x-coordinates of flood delta polygon
    S.yfloodpol=[];                                                            % y-coordinates of flood delta polygon
    S.ldbspit='';                                                              % file with x and y coordinates of the spit polygon [Nx2] or directly provide a matrix
    S.xspitpol=[];                                                             % x-coordinates of the spit polygon
    S.yspitpol=[];                                                             % y-coordinates of the spit polygon
    S.dxf=50;                                                                  % resolution of flood delta area [m]
    S.overdepth=2;                                                             % initial overdepth flood delta [m]
    %% ------------------------ formatting / plotting -------------------------
    S.plotvisible=1;                                                           % plot and update figure with wave conditions and modeled shoreline during run
    S.xlimits=[];                                                              % x-limits of plot [m] [1x2] <- leave empty to automatically do this
    S.ylimits=[];                                                              % y-limits of plot [m] [1x2] <- leave empty to automatically do this
    S.xywave =[];                                                              % X,Y location and scale of wave arrow [m] [1x2] (automatically determined on the basis of xy-limits and wave angle
    S.xyoffset=[0,0];                                                          % shift in X,Y locaton for plotting <- leave empty to automatically shift grid <- use [0,0] for no shift
    S.pauselength=[];                                                          % pause between subsequent timesteps (e.g. 0.0001) <- leave empty to not pause plot
    S.ldbplot = {};                                                            % cell array with at every line a string with LDB-filename, string with legend entry, string with plot format (e.g. 'b--') <- e.g. {'abc.ldb','line 1','k--'; 'def.ldb','line 2','r-.'; etc} <- leave empty to not use additional plots
    S.ploths = 0;                                                              % plot wave height at depth-of-closure (TDP) as coloured markers and text along the coast (use 0/1 as switch). A larger value than 1 plots at every 'nth' grid cell (so 5 at every five cells).
    S.plotdir = 0;                                                             % plot wave direction at depth-of-closure (TDP) as quivers and text along the coast (use 0/1 as switch). A larger value than 1 plots at every 'nth' grid cell (so 5 at every five cells).
    S.plotqs = 0;                                                              % plot transport rates at depth-of-closure (TDP) as coloured markers and text along the coast (use 0/1 as switch). A larger value than 1 plots at every 'nth' grid cell (so 5 at every five cells).
    S.plotupw = 0;                                                             % plot locations with high angle correction (use 0/1 as switch)
    S.llocation='SouthWest';                                                   % location of legend. Shortened for Octave compatibility
    S.ld=3000;                                                                 % Width of the land fill behind the shoreline [m]
    S.usefill = 1;                                                             % option switch that can be used to only plot lines instead of the fill
    S.usefillpoints = 0;                                                       % option that can be used to force the model to use a specified number of landpoints landward of open coastlines for the plotting of filled-land (e.g. 6 means that a landward points is computed for 6 points, while 0 means the model automatically determines a land fill at 1 of the 4 sides)
    S.fignryear=12;                                                            % the number of plots that need to be stored to a file each year [1/year]
    S.plotinterval=1;                                                          % the interval at which the coastline is plotted (in number of timesteps inbetween)
    S.figplotfreq=[];                                                          % stores the coastline at this frequency w.r.t. model start, and plots it [years of interval]
    S.fastplot=1;                                                              % use imwrite plotter, that is ~2x as fast as print
    %% -------------------------------- output --------------------------------
    S.outputdir='Output\';                                                     % output directory
    S.outputfile = 'shorelines_output';                                        % default filename (either the .mat or .nc extension will be added)
    S.rundir='Delft3D\def_model\';                                             % run directory of a coupled Delft3D flow model
    S.xyout=[];                                                                % specify multiple curvi-linear output grids that remain fixed at storageinterval frequency. There are 2 options: Option 1: Specify two coordinates as [x1,y1,x2,y2]. A linear grid is made with ds0 spacing. Each line in the matrix represents a new grid. Option 2: Specify a grid as a cell {} with [Nx2] arrays of xy-coordinates. The grid is used 1 on 1. You can specify multiple grids by adding more cells. For example: {[5x2],[12x2],...} 
    S.xyprofiles=[];                                                           % specify multiple output profile locations [Nx2 with x,y; ...]. The model will determine the closest coastline and dune point. 
    S.storageinterval=50;                                                      % time interval of storage of output data to a file ('output.mat'; in [day])
    S.storagedate={};                                                          % moments in time at which data needs to be stored {'yyyy-mm-dd', 'yyyy-mm-dd', 'yyyy-mm-dd', ...}, when specified this keyword is used instead of 'storageinterval'
    S.netcdf=0;                                                                % switch for using nc-files (1) or the mat-files (0) as output
    S.separatepgrids=1;                                                        % creates separate nc-files for each of the P-grids (when 1) while the p-grids will be included as a group in the main ncfile otherwise (when 0)
    %% -------------------------- extract shorelines --------------------------
    S.slplot={};                                                               % slplot
    S.extractxy=0;                                                             % get file with shorelines coordinates x,y
    S.printfig=0;                                                              % switch to print figures to files (0/1)
    %% ---------------- extract shoreline & dune foot locations ---------------
    S.yesplot=0;                                                               % only needed for interactive mode, not for regular runs with plots
    S.bermwplot=0;                                                             % for extracting  against certain times for a certain transect(s).
    S.bermwplotint=[];                                                         % beach berm width plot interval [Months]
    S.qplot=0;                                                                 % to plot wave and wind transport at each time step for a certain transect(s).
    S.transect='';                                                             % file indictes x-y transects to be plotted
    S.clplot=0;                                                                % to track coastline location relative to the initial coasltime against certain time interval for a certain transect(s).
    S.clplotint=[];                                                            % coastline change plot interval a certain transect(s) [Months] .
    S.extractbermplot=0;                                                       % to allow data extraction for berm width plotting
    %% video 
    S.video=0;                                                                 % switch for capturing a video (0/1)
    %% debug
    S.debug=0;                                                                 % switch for setting the debug modus (0=off, 1-5 different storage / plot options)
    %% --------------------------- data Assimilation---------------------------
    S.da=0;                                                                    % data assimilation switch (0/1)
    S.bs=0;	                                                                   % switch of brier-skill score (0/1)

    %% To handle Octave runs with incompatible mex function, set global handle to compatible intersection function
    
    %% SUBSITUTE MODEL INPUT IN THE S STRUCTURE, AND USE DEFAULTS IF NO VALUE IS SPECIFIED
    %if (isoctave)   
    %   get_intersections=@get_intersections;
    %else
    %   get_intersections=@mexinterx;
    %end
    fieldnms=get_fields(S0);
    for ii=1:length(fieldnms)
        % make sure that the 'underscores' are removed from input keywords, as backward compatability option
        if ~isempty(findstr(fieldnms{ii},'_')) || sum(lower(fieldnms{ii})~=fieldnms{ii})~=0
            fieldnew=regexprep(fieldnms{ii},'_','');
            fieldnew=lower(fieldnew);
            S0.(fieldnew)=S0.(fieldnms{ii});
            S0=rmfield(S0,fieldnms{ii});
            fieldnms{ii}=fieldnew;
        end
    
        % check if fieldname exists, otherwise provide a warning.
        if ~isfield(S,fieldnms{ii}) && ~strcmpi(fieldnms{ii},'x_mc') && ~strcmpi(fieldnms{ii},'y_mc')
            fprintf('   - Warning : ''%s'' is not a default keyword.\n',fieldnms{ii});
        end
        S.(fieldnms{ii}) = S0.(fieldnms{ii});
    end
end

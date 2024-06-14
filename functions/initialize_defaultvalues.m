function [S]=initialize_defaultvalues(S0)
% function [S]=initialize_defaultvalues(S0)
%
% INPUT
%    S       Default settings of the ShorelineS model (as specified in this function)
%    S0      User input given to the ShorelineS model
%
% OUTPUT
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
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses
%   --------------------------------------------------------------------

    fprintf('%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');  
    fprintf('%s\n','%                  RUN SHORELINES                   %');  
    fprintf('%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');  
    fprintf('%s\n','  Initialize default values');

    if (isoctave)
        warning('off','all');
    end
    
    %% ---------------------- simulation wave parameters ----------------------        
    S.Hso=1;                                                                   % wave height [m]
    S.tper=6;                                                                  % peak wave period [s]
    S.phiw0=330;                                                               % deep water wave angle in degrees [?N]
    S.spread=90;                                                               % wave spreading [?] (wave_dir from range:  S.phiw0 +/- 0.5*S.spread)
    S.WVCfile='';                                                              % wave time-series <-leave empty to use wave parameters ('S.Hso', 'S.phiw0' and 'S.spread')
    S.ddeep=25;                                                                % Waterdepth the location of wave climate, corresponding with S.Hso [m]
    S.dnearshore=8;                                                            % Waterdepth at the 'dynamic boundary', corresponding with S.phif [m]
    S.randomseed=-1;                                                           % Seed to generate randomized series
    S.interpolationmethod='weighted_distance';  %'alongshore_mapping';         % Method for interpolating wave data of multiple wave stations on the coast (either 'alongshore_mapping' / 'weighted_distance')
    %% ------------------- simulation coastline definition --------------------        
    S.LDBcoastline='';                                                         % LDB with initial coastline shape ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
    S.gis_convention = 0;                                                      % coastline defined accordsing to GIS convention (e.g. counter clockwise)
    S.ds0=100;                                                                 % initial space step [m]
    S.griddingmethod=2;                                                        % method for regenerating the grid (1: only splitting and merging cells if they are too small, 2: uniform grid regeneration if criteria for gridsize or exceeded)
    S.d=10;                                                                    % active profile height [m]
    S.phif=[];                                                                 % Orientation of the foreshore [?N] (in degrees) <- only relevant for 'KAMP', 'MILH' or 'VR14'
    S.maxangle=60;                                                             % maximum coastline re-orientation between individual grid cells (affecting spit width and stabilizing small scale features in case of dense grids)
    %% ------------------- simulation transport parameters --------------------    
    S.trform='CERC';                                                           % switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
    S.b=1e6;                                                                   % CERC : coeff in simple cerc formula
    S.qscal=1;                                                                 % Calibration factor of the transport (works for all transport formulas)
    S.d50=2.0e-4;                                                              % KAMP & MILH & VR14 : median grain diameter [m]
    S.d90=3.0e-4;                                                              % KAMP & MILH & VR14 : median grain diameter [m]
    S.porosity=0.4;                                                            % KAMP & MILH & VR14 : S.porosity (typically 0.4) [-]
    S.tanbeta=0.03;                                                            % KAMP & MILH & VR14 : mean bed slope [ratio 1/slope]
    S.rhos=2650;                                                               % KAMP & MILH & VR14 : density of sand [kg/m3]
    S.rhow=1025;                                                               % KAMP & MILH & VR14 : density of water [kg/m3]
    S.g=9.81;                                                                  % KAMP & MILH & VR14 : gravitational acceleration [m2/s]
    S.alpha=1.8;                                                               % KAMP & MILH & VR14 : calibration factor for point of breaking (S.alpha = 1.8 for Egmond data)
    S.gamma=0.72;                                                              % KAMP & MILH & VR14 : breaking coefficient (Hs/h) with 5% breaking waves
    S.Pswell=20;                                                               % VR14 : Percentage swell (between 0 - 100) [-]
    S.crit=.9;                                                                 % stability criterion (not active)
    S.ks=0.05;                                                                 % 
    S.hclosure=8;                                                              % 
    S.Cf=0.0023;                                                               % roughness factor
    S.hmin=0.1;                                                                % 
    S.sphimax=[];                                                              % the computation of the maximum angle of the waves can be bypassed by specifying this sphimax value (e.g. at sphimax=42) or by using 'auto' which means that it is computed only at t0
    S.relaxationlength=[];                                                     % length over which transport decelerates in meters, which adds inertia to the longshore current. It scales linearly with the wave height below 1m waves.
    S.suppress_highangle=0;                                                    % switch 0/1 to disable the high-angle instabilities by limiting the transport angle to the critical high-angle orientation (when it is set at 1)
    %% ------------------- simulation time steps & numerical-------------------
    S.tc=1;                                                                    % switch for using adaptive time step. The tc-value gives a fraction of the automatically computed timestep. So when S.tc=0.5 is used, it will use a timestep which is half of the maximum step that is permissable. S.tc=0 means that the fixed step S.dt is used instead.
    S.dt=0;                                                                    % value of a fixed timestep for the coastline processes [fraction of a year]. It is only used when S.tc=0;
    S.dtdune=[];                                                               % value of a fixed timestep for dune process [fraction of a year]. It should always be less than 'dt' of the coastline;
    S.reftime='2020-01-01';                                                    % reference/start time of the simulation ('yyyy-mm-dd') <- leave empty to use t=0
    S.endofsimulation='2040-01-01';                                            % the end time of the simulation ('yyyy-mm-dd')
    S.twopoints=1;                                                             % switch for 'S.twopoints approach' which determines the type of response during high-angle wave events (default=1)
    S.smoothfac=0;                                                             % smoothing factor used to re-arrange grid every timestep (only for griddingmethod==1), with reasonable values between 0 and 0.1.
    S.smoothrefrac=0;                                                          % smoothing fraction of the coastline orientation (PHIcs) which is used only for the refraction of waves from nearshore location (TDP) to point of breaking (BR), and not for the transport computation. Values between 0 and 1.
    %% -------------------------- boundary condition -------------------------- for (non cyclic) sections (ex.straight shoreline) ,If more than one (non cyclic)sections should adjusted manually
    S.boundary_condition_start='Fixed';                                        % left boundary condition (either 'Fixed', 'Closed', 'Angleconstant' or 'Periodic'). An additional argument with the specified coastangle or transport rate can be given by providing a cell-format. For example, {'Closed',-6000} to specify a negative transport of -6000 m3/yr or {'Angleconstant',300} to enforce a coastangle of 300캮. 
    S.boundary_condition_end='Fixed';                                          % right boundary condition (either 'Fixed', 'Closed', 'Angleconstant' or 'Periodic'). An additional argument with the specified coastangle or transport rate can be given by providing a cell-format. For example, {'Closed',-6000} to specify a negative transport of -6000 m3/yr or {'Angleconstant',300} to enforce a coastangle of 300캮. 
    %% ------------------------ climate change impact -------------------------
    S.ccSLR=[];                                                                % climate impacted rise in sea level (SLR). This can be a constant rate of sea level rise (e.g. S.ccSLR=0.002 m/yr) or a table with the absolute sea level against time [Nx2]. Wherein 'time in datenum format' and 'sea level with respect to initial situation'). Note that the rates per year are computed automatically! For example, the rate of SLR goes up from 2mm/yr to 10mm/yr in this table: S.ccSLR=[t1, 0.002; t2,0.003; t3,0.005; t4,0.01]; The S.tanbeta is used as 'slope angle' for the BRUUN rule. 
    S.ccHS=[];                                                                 % climate impacted increase of the wave height (HS). This can be a constant rate per year (e.g. S.ccHS=0.001, which is +0.1% increase in HS per year) or a table with for various time instances the relative change in wave height with respect to the initial situation [Nx2]. Wherein 'time in datenum format' and 'relative change in wave height w.r.t. initial situation' as a factor. For example, the HS may increase by up to 40% up till 't4' in this table S.ccHS=[t1, 0.0; t2,0.1; t3,0.25; t4,0.4]; 
    S.ccDIR=[];                                                                % climate impacted change in wave direction (DIR). This can be a constant rate (e.g. S.ccDIR=0.05 째/yr) or a table with the wave direction change with respect to the initial situation against time [Nx2]. Wherein 'time in datenum format' and 'relative change in wave direction w.r.t. initial situation' as a # degrees. For example, the wave direction may be adjsuted by up to 0째 up till +3째 at 't4' in this table S.ccDIR=[t1, 0째; t2,+0.5째; t3,+1.5째; t4,+3째]; 
    %% ----------------------------- structures -------------------------------
    S.struct=1;                                                                % switch for using hard structures
    S.LDBstructures='';                                                        % LDB with hard structures ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
    S.x_hard=[];                                                               % x-coordinates of hard structures ([1xN]) seperated by nans <- leave empty to use interactive mode!
    S.y_hard=[];                                                               % y-coordinates of hard structures ([1xN]) seperated by nans <- leave empty to use interactive mode!
    S.structtype={};
    S.Aw=1.27;                                                                 % factor for determining depth of closure at bypassing groyne (1.27 if time series is used) This value is used by default.
    S.Awfixedhs=5;                                                             % factor for determining depth of closure at bypassing groyne ir a representative Hs is used instead of a climate or timeseries. This value is used instead of 'Aw' if S.WVCfile is empty.
    S.bypasscontractionfactor=1;                                               % scaling factor for the transport bypass at a groyne as a result of contraction of the flow (always >=1). Setting the bypass contraction factor larger than 1 means that the accretion does not go to the tip of the structure. 
    S.bypassdistribution_power=1;                                              % this power controls the distribution of bypassed sediment of groynes in the shadow zone. It will be a triangle with pwr=1 (default), and more sediment will end up closer to the structure for a higher power (e.g. pwr=2).
    %% ------------- wave transmission over submerged breakwater --------------
    S.transmission=0;                                                          % Switch to allow for wave transmission over breakwaters
    S.transmform='angr';                                                       % Wave transmission approaches for rough structures: d'Angremond (1996) 'angr', Van Gent et al. for rough structures (2023) 'gent' or Seabrook and Hall (1998) 'seabrhall'        
    S.transmdir=1;                                                             % Switch to adapt wave direction after transmission over breakwater scaling with percentage transmission [S.transmdir=1], or keep original wave direction [S.transmdir=0];
    S.transmfile='';                                                           % File with characteristics of submergable breakwaters ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
    S.transmbwdepth=[];                                                        % depth at which breakwater is positioned [m] Positive downward.
    S.transmcrestheight=[];                                                    % breakwater height [m] w.r.t. MSL Positive upward, above MSL.
    S.transmslope=[];                                                          % breakwater slope [-] 
    S.transmcrestwidth=[];                                                     % breakwater crest width [m] 
    S.transmd50=1;                                                             % D50 of armour layer, needed as input for Seabrook and Hall formulation 'SeabrHall'. Put to 1 m as default.                
    %% ------------------------- permeable structures -------------------------
    S.perm=0;                                                                  % switch for using hard structures
    S.LDBpermeable='';                                                         % LDB with perm structures ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
    S.x_perm=[];                                                               % x-coordinates of permeable structures ([1xN]) seperated by nans <- leave empty to use interactive mode!
    S.y_perm=[];                                                               % y-coordinates of permeable structures ([1xN]) seperated by nans <- leave empty to use interactive mode!
    S.wavetransm=[1]; 
    %% ----------------------------- revetments -------------------------------
    S.revet=1;                                                                 % switch for using revetments
    S.LDBrevetments='';                                                        % LDB with hard structures ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
    S.x_revet=[];                                                              % x-coordinates of revetments ([1xN]) seperated by nans <- leave empty to use interactive mode!
    S.y_revet=[];                                                              % y-coordinates of revetments ([1xN]) seperated by nans <- leave empty to use interactive mode!
    S.iterrev=5;                                                               % number of iterations used to determine bypassing along a revetment (sometimes with very long revetments, with spatially varying conditions a large value is needed for a stable computation)
    S.crit_width=5;                                                            % Critical width used for determing coastal width in front of revetments. This scales the bypass accordingly!
    %% --------------------------- wave diffraction ---------------------------
    S.diffraction=0;                                                           % wave diffraction < use 1 to avtivate
    S.wdform='Roelvink';                                                       % Wave diffraction approach to treat angles (Roelvink(default), Hurst). With Hurst the omegat factor is set to 0 and the rotfac to 1, and directional spreading is handled with discrete directional sectors later in the wave_difffraction_coeff function. 
    S.kdform='Kamphuis';                                                       % Computation of kd according to 'Kamphuis' or 'Roelvink' analytical approx. The Roelvink method has been derived on the basis of our data. The Kamphuis method gives a smoother result for the wave energy distribution. 
    S.rotfac=1.5;
    S.diffdist=[];                                                             % Specify the maximum distance of structures to coastline points with diffraction (limits influence area)
    %% ---------------------------- nourishments ------------------------------
    S.nourish=0;                                                               % switch (0/1) for nourishments
    S.growth=1;                                                                % calibration of nourishment growth rate
    S.norfile='';                                                              % file with nourishments, which stores a table with at each line the properties of a nourishment (x1,y1,x2,y2,t1,t2,volume)
    S.LDBnourish='';                                                           % LDB with nourishment locations ([Nx2] ASCII FILE WITHOUT HEADER) (i.e. polygon around relevant grid cells) <- leave empty to use interactive mode!
    S.nourratefile='';                                                         % nourishment rates placed in order for each nourishment polygons
    S.nourstartfile='';                                                        % nourishment start dates placed in order for each  nourishment polygons
    S.nourendfile='';                                                          % nourishment end dates placed in order for each  nourishment polygons
    S.nourrate=100;                                                            % 
    S.nourstart=0;                                                             % 
    S.nourend=[];                                                              % 
    %% ------------------------ shoreface nourishments ------------------------
    S.fnourish=0;                                                              % switch (0/1) for shoreface nourishments
    S.fnorfile='';                                                             % shoreface nourishment input file 
    S.mb=-0.5;                                                                 % 
    S.labda0=0.56E-6;                                                          % 
    S.rhos=2650;                                                               % 
    S.sal=35;                                                                  % 
    S.temp=10;                                                                 % 
    S.K=[];                                                                    % 
    %% -------------------------- sources and Sinks ---------------------------
    S.sources_sinks='';                                                        % 
    S.SSfile='';                                                               % 
    %% -------------------- aeolian transport to the dunes --------------------
    S.dune=0;                                                                  % switch for using estimate dune evolution
    S.LDBdune='';                                                              % .dun filename of file with dune parameters [Nx5] with xdune,ydune,wberm,Dfelev,Dcelev in the columns (ASCII FILE WITHOUT HEADER)  <- leave empty to use default values 
    S.kf=0.02;                                                                 % Friction coefficient
    S.Cs=5e-4;                                                                 % Erosion coefficient of the dunes during storms, which scales the rate of erosion
    S.Cstill=5e-6;                                                             % Erosion coefficient of the dunes during storms for dunes with consolidated till layers, which scales the rate of erosion
    S.xtill=[];                                                                % Thickness of layer with sand at the seaward side of the dunes (which can be eroded first) before a till layer is exposed. Leave empty (default) if dunes are 100% sandy (i.e. no till layers), or use a very large thickness then. A spatially varying value can be specified by using a [Nx3] table with x,y,xtill per location at each line. 
    S.perctill=80;                                                             % Percentage of till that is present in the till layer of the dunes (i.e. when dune erosion is beyond xtill). A spatially varying value can be specified by using a [Nx3] table with x,y,perctill per location at each line. 
    S.d50r=2.5e-4;                                                             % Median reference grain size
    S.rhoa=1.225;                                                              % Air density []
    S.duneAw=0.1;                                                              % Coefficient (Bagnold, 1937)
    S.Kw=4.2;                                                                  % Empirical coefficient (Sherman et al. 2013)
    S.k=0.41;                                                                  % Von Karman's coefficient
    S.segmaw=0.1;                                                              % Empirical factor used for scaling impact of the fetch length
    S.maxslope=1/15;                                                           % The maximum slope angle (1:slope) with default of 1:20. The dunefoot height is lowered if the beach gets too steep (preserving the max slope).
    S.A_overwash=3;
    S.xdune=0;                                                                 % x coords where Wberm and Dfelev are given (only used if LDBdune is empty)
    S.ydune=0;                                                                 % y coords where Wberm and Dfelev are given (only used if LDBdune is empty)
    S.Wberm=50;                                                                % Berm width (m), relevant for the fetch of the eolian transport & the run-up of Stockdon and Ghonim (only used if LDBdune is empty)
    S.Dfelev=3;                                                                % Dune foot elevation w.r.t. MSL (m) (only used if LDBdune is empty)
    S.Dcelev=8;                                                                % Dune crest elevation w.r.t. MSL (m) (only used if LDBdune is empty)
    %% --------------------------- wind conditions ----------------------------
    S.WNDfile='';                                                              % wind time-series filename <-leave empty to use 'S.uz', 'S.phiwnd0' and 'S.spread'
    S.Cd=0.002;                                                                % wind drag coefficient (-)
    S.uz='';                                                                   % wind velocity at z (m)
    S.z=10;                                                                    % elevation of measured wind data
    S.phiwnd0=330;                                                             % wind angle [degN]
    %% --------------------- still water levels and run-up---------------------
    S.runupform='Stockdon';                                                    % switch for run-up formulation (default 'Stockdon', or 'Ghonim', or 'Larson')
    S.runupfactor=1;                                                           % tuning factor for runup (linear)
    S.WATfile='';                                                              % Water levels time series file with [Nx2] date/time in 'yyyymmddHHMM' and waterlevel relative to MSL
    S.WVDfile='';                                                              % Wave height time series file with [Nx2] date/time in 'yyyymmddHHMM' and wave height, wave period and wave direction
    S.SWL0=0;                                                                  % Fixed still water level relative to MSL
    %%------------------------- Sediment limitations --------------------------
    S.sedlim=1;                                                                % Switch for sediment limiter
    S.LDBsedlim='';                                                            % file with sediment limitation coordinates and properties, with [Nx3] specification. The columns have the x and y coordinates and width at which transport starts to reduce gradually to zero, separated with NaN's for different sections.
    S.x_sedlim=[];                                                             % alternative input with specifying x-coordinates of areas with sediment limitation separately, as [Nx1] separated with NaN's for different sections
    S.y_sedlim=[];                                                             % alternative input with specifying y-coordinates of areas with sediment limitation separately, as [Nx1] separated with NaN's for different sections
    S.width_sedlim=[];                                                         % alternative input with specifying width at which transport starts to reduce for areas with sediment limitation separately, as [Nx1] separated with NaN's for different sections
    %% --------------------------- mud properties -----------------------------
    S.mud=0;                                                                   % mud transport option
    S.mud_taucr=0.3;                                                           % critical shear stress for erosion (N/m2)
    S.mud_M=1.e-4;                                                             % erosion rate (kg/m2/s)
    S.mud_B=1000;                                                              % muddy transport zone width (m)
    S.mud_MHW=1;                                                               %
    S.mud_Tfm=10;                                                              %
    S.cyclic=0;                                                                % switch for cyclic boundary conditions
    %% ------------------- physics of spit width and channel ------------------
    S.spit_method='default';                                                   % Overwash method
    S.spit_width=50;                                                           % width of tip of spit (used for overwash)
    S.spit_headwidth=200;                                                      % width of tip of spit (used for upwind correction)
    S.OWscale=0.1;                                                             % scales the rate of the overwash per timestep (i.e. what part of the deficit is moved to the backbarrier)
    S.OWtimescale=0.0;                                                         % timescale for overwash (i.e. what part of the deficit is moved to the backbarrier)
    S.spit_Dsf=S.d*0.8;                                                        % underwater part of active height for shoreface -> used only in spit-width function
    S.spit_Dbb=0.5*S.spit_Dsf;                                                 % underwater part of active height for back-barrier -> used only in spit-width function
    S.Bheight=2;                                                               % berm height used for overwash funciton (i.e. added to Dsf or Dbb)
    S.tide_interaction=false;                                                  % 
    S.wave_interaction=false;                                                  % 
    S.wavefile='';                                                             % wave table (.mat)
    S.surf_width_w=250;                                                        % width of surf zone, where to collect the wave conditions from wave table
    S.surf_width=1000;                                                         % width of surf zone, where to update the bathymetry
    S.bathy_update='';                                                         % the dates when the bathymetry should be updated, the input should be in dates form, can accept more than one  {'yyyy-mm-dd'};
    %% ------------------------------- channel --------------------------------
    S.channel=0;                                                               % switch (0/1)for migrating inlet on
    S.channel_width=550;                                                       % target channel width
    S.channel_fac=0.08;                                                        % adaptation factor
    S.channel_disch_rate=0;                                                    % adaptation factor
    S.channel_disch_R=300;                                                     % adaptation factor
    S.LDBchannel=[];                                                           % initial channel axis
    S.flood_delta=0;                                                           % switch (0/1) for flood delta losses
    S.LDBflood=[];                                                             % wide outline of potential flood delta deposits
    S.x_flood_pol=[];                                                          % 
    S.y_flood_pol=[];                                                          % 
    S.dxf=50;                                                                  % resolution of flood delta area
    S.overdepth=2;                                                             % initial overdepth flood delta
    S.R_ero=300;                                                               % radius of flood delta influence on coast
    S.R_depo=600;                                                              % radius of inlet influence on flood delta
    S.Tscale=1;                                                                % timescale of flood delta filling
    S.xr_mc='';                                                                % 
    S.yr_mc='';                                                                % 
    %% ------------------------- formatting / output --------------------------
    S.plotvisible=1;                                                           % plot and update figure with wave conditions and modeled shoreline during run
    S.xlimits=[];                                                              % X limits of plot [m] [1x2] <- leave empty to automatically do this
    S.ylimits=[];                                                              % Y limits of plot [m] [1x2] <- leave empty to automatically do this
    S.XYwave =[];                                                              % X,Y location and scale of wave arrow [m] [1x2] (automatically determined on the basis of xy-limits and wave angle
    S.XYoffset=[0,0];                                                          % shift in X,Y locaton for plotting <- leave empty to automatically shift grid <- use [0,0] for no shift
    S.pauselength=[];                                                          % pause between subsequent timesteps (e.g. 0.0001) <- leave empty to not pause plot
    S.outputdir='Output\';                                                     % output directory
    S.rundir='Delft3D\def_model\';                                             % run directory of a coupled Delft3D flow model
    S.LDBplot = {};                                                            % cell array with at every line : string with LDB-filename, string with legend entry, string with plot format (e.g. 'b--') <- e.g. {'abc.ldb','line 1','k--'; 'def.ldb','line 2','r-.'; etc} <- leave empty to not use additional plots
    S.plotHS = 0;                                                              % plot wave height at depth-of-closure (TDP) as coloured markers and text along the coast (use 0/1 as switch). A larger value than 1 plots at every 'nth' grid cell (so 5 at every five cells).
    S.plotDIR = 0;                                                             % plot wave direction at depth-of-closure (TDP) as quivers and text along the coast (use 0/1 as switch). A larger value than 1 plots at every 'nth' grid cell (so 5 at every five cells).
    S.plotQS = 0;                                                              % plot transport rates at depth-of-closure (TDP) as coloured markers and text along the coast (use 0/1 as switch). A larger value than 1 plots at every 'nth' grid cell (so 5 at every five cells).
    S.plotUPW = 0;                                                             % plot locations with high angle correction (use 0/1 as switch)
    S.llocation='SouthWest';                                                   % location of legend. Shortened for Octave compatibility
    S.ld=3000;                                                                 % Width of the land fill behind the shoreline [m]
    S.usefill = 1;                                                             % option switch that can be used to only plot lines instead of the fill
    S.usefillpoints = 0;                                                       % option that can be used to force the model to use a specified number of landpoints landward of open coastlines for the plotting of filled-land (e.g. 6 means that a landward points is computed for 6 points, while 0 means the model automatically determines a land fill at 1 of the 4 sides)
    S.fignryear=12;                                                            % the number of plots that need to be stored to a file each year [1/year]
    S.plotinterval=1;                                                          % the interval at which the coastline is plotted (in number of timesteps inbetween)
    S.fastplot=1;                                                              % Use imwrite plotter, that is ~2x as fast as print
    S.xyout=[];                                                                % Specify multiple curvi-linear output grids that remain fixed at storageinterval frequency. There are 2 options: Option 1: Specify two coordinates as [x1,y1,x2,y2]. A linear grid is made with ds0 spacing. Each line in the matrix represents a new grid. Option 2: Specify a grid as a cell {} with [Nx2] arrays of xy-coordinates. The grid is used 1 on 1. You can specify multiple grids by adding more cells. For example: {[5x2],[12x2],...} 
    S.xyprofiles=[];                                                           % Specify multiple output profile locations [Nx2 with x,y; ...]. The model will determine the closest coastline and dune point. 
    S.storageinterval=50;                                                      % Time interval of storage of output data to a file ('output.mat'; in [day])
    %% -------------------------- extract shorelines --------------------------
    S.SLplot={};                                                               % 
    S.extract_x_y=0;                                                           % get file with shorelines coordinates x,y
    S.print_fig=0;                                                             % 
    %% ---------------- extract shoreline & dune foot locations ---------------
    S.yesplot=0;                                                               % only needed for interactive mode, not for regular runs with plots
    S.bermw_plot=0;                                                            % For extracting  against certain times for a certain transect(s).
    S.bermw_plot_int=[];                                                       % beach berm width plot interval [Months]
    S.qplot=0;                                                                 % to plot wave and wind transport at each time step for a certain transect(s).
    S.transect='';                                                             % file indictes x-y transects to be plotted
    S.CLplot=0;                                                                % to track coastline location relative to the initial coasltime against certain time interval for a certain transect(s).
    S.CLplot_int=[];                                                           % coastline change plot interval a certain transect(s) [Months] .
    S.extract_berm_Plot=0;                                                     % to allow data extraction for berm width plotting
    %% video
    S.video=0;                                                                 % 
    %% debug
    S.debug=0;                                                                 % 
    %% --------------------------- data Assimilation---------------------------
    S.DA=0;                                                                    % Data assimilation Switch
    S.BS=0;	

    %% To handle Octave runs with incompatible mex function, set global handle to compatible intersection function
    
    %% SUBSITUTE MODEL INPUT IN THE S STRUCTURE, AND USE DEFAULTS IF NO VALUE IS SPECIFIED
    %if (isoctave)   
    %   get_intersections=@get_intersections;
    %else
    %   get_intersections=@mexinterx;
    %end
    fieldnms=get_fields(S0);
    for ii=1:length(fieldnms)
        if ~isfield(S,fieldnms{ii}) && ~strcmpi(fieldnms{ii},'x_mc') && ~strcmpi(fieldnms{ii},'y_mc')
            fprintf('   - Warning : ''%s'' is not a default keyword.\n',fieldnms{ii});
        else
            S.(fieldnms{ii}) = S0.(fieldnms{ii});
        end
    end
end

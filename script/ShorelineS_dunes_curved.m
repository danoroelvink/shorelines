addpath('..\functions\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;
S.twopoints=1; 
S.Hso=1.5;                                                                 % wave height [m]
S.phiw0=270;                                                               % deep water wave angle [°N]
S.spread=0;                                                                % wave spreading [°] (wave_dir from range:  S.phiw0 +/- 0.5*S.spread)
S.d=8;                                                                     % active profile height [m]
S.ddeep=56;                                                                % water depth at deep water (For refraction)
S.dnearshore=4;                                                            % water depth at nearshore (For refraction)
S.LDBcoastline='data\curved_coastline.ldb';                                % LDB with initial coastline shape ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
S.trform='CERC3';                                                          % switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
S.b=2.9833E+6;                                                             % CERC : coeff in simple cerc formula
S.tper=9.456;                                                              % KAMP & MILH : peak wave period [s]
S.d50=2.5e-4;                                                              % KAMP & MILH & VR14 : median grain diameter [m]
S.porosity=0.4;                                                            % KAMP & MILH & VR14 : S.porosity (typically 0.4) [-]
S.tanbeta=0.05;                                                            % KAMP & MILH & VR14 : mean bed slope [ratio 1/slope]
S.rhos=2650;                                                               % KAMP & MILH & VR14 : density of sand [kg/m3]
S.rhow=1025;                                                               % KAMP & MILH & VR14 : density of water [kg/m3]
S.g=9.81;                                                                  % KAMP & MILH & VR14 : gravitational acceleration [m2/s]
S.alpha=1.8;                                                               % KAMP & MILH & VR14 : calibration factor for point of breaking (S.alpha = 1.8 for Egmond data)
S.gamma=0.8;                                                               % KAMP & MILH & VR14 : breaking coefficient (Hs/h) with 5% breaking waves
S.Pswell=80;                                                               % VR14 : Percentage swell (between 0 - 100) [-]
S.crit=.9;                                                                 % stability criterion (not active)
S.reftime='2020-01-01';                                                    % Reference time (i.e. 'yyyy-mm-dd') <- leave empty to use t=0
S.endofsimulation='2020-03-01';                                            % time step [year] (use 1/365/4 for time-series) (not used in combination with adaptive time step)
S.ds0=100;                                                                 % initial space step [m]
S.tc=0;
S.dt=1/365/2;                                                              % 12 hours 
S.smoothfac=0;
S.Courant=1;
S.ns=50;                                                                   % number of ... (not used)
S.struct=0;                                                                % switch for using hard structures
S.LDBstructures='';                                                        % LDB with hard structures ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
S.growth=0;                                                                % calibration of nourishment growth rate
S.nourish=0;                                                               % switch (0/1) for nourishments
S.LDBnourish='';                                                           % LDB with nourishment locations ([Nx2] ASCII FILE WITHOUT HEADER) (i.e. polygon around relevant grid cells) <- leave empty to use interactive mode!
S.nourrate=100;
S.nourstart=0;
S.OWscale=0.2;
S.xlimits=[-100,1200];                                                     % X limits of plot [m] [1x2] <- leave empty to automatically do this                                                          % Y limits of plot [m] [1x2] <- leave empty to automatically do this
S.ylimits=[-1600,1600];
S.XYwave =[];                                                              % X,Y location and scale of wave arrow [m] [1x2] (automatically determined on the basis of xy-limits and wave angle
S.XYoffset=[0,0];                                                          % shift in X,Y locaton for plotting <- leave empty to automatically shift grid <- use [0,0] for no shift
S.pauselength=[];                                                          % pause between subsequent timesteps (e.g. 0.0001) <- leave empty to not pause plot
S.fignryear=5;
S.storageinterval=30;                                                      % Time interval of storage of output file ('output.mat'; [day])
S.outputdir='Output\dunes_curved\';                                        % output directory for plots and animation
S.usefill = 1;
S.usefillpoints = 6;
S.boundary_condition_start='Closed';                                       % boundary condition: 'Angleconstant', 'Fixed' , 'Closed', 'Periodic' , see "transport_boundary_condition.m"
S.boundary_condition_end='Closed';
S.randomseed = 100;

%% Extra inputs for dune evolution simulation
S.dune = 1; 
S.RWSfiletype=true;
S.LDBdune='data\curved_dunefoot_3m.dun';                                   % LDB with initial dune foot position,berm width, dune foot height and dune crest height   [Nx2] 
S.uz=15;                                                                   % wind velocity at z (m)
S.z=10;                                                                    % elevation of measred wind data
S.phiwnd0=270;                                                             % wind angle [°N]
S.SWL0=1.5;                                                                % Fixed still water level relative to MSL
S.Cs=7.5e-5;                                                               % Impact coefficient waves based on thesis M. Ghonim || very sensitive
S.runupform='Stockdon'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN SHORELINES MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
[S,O]=ShorelineS(S);

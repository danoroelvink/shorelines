addpath(genpath('..\functions\'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;
S.smoothfac=0;
S.Hso=1;                                                                   % wave height [m]
S.phiw0=235;                                                               % deep water wave angle [°N]
S.spread=0;                                                                % wave spreading [°] (wave_dir from range:  S.phiw0 +/- 0.5*S.spread)
S.d=5;                                                                     % active profile height [m]
S.ddeep=100;                                                               % water depth at deep water (For refraction)
S.dnearshore=100;                                                          % water depth at nearshore (For refraction)
S.LDBcoastline='Data\spit.xy';                                             % LDB with initial coastline shape ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
S.phif=270;                                                                % Orientation of the foreshore [°N] <- only relevant for 'KAMP', 'MILH' or 'VR14'
S.trform='VR14';                                                           % switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
S.b=1E+6;                                                                  % CERC : coeff in simple cerc formula
S.tper=7;                                                                  % KAMP & MILH : peak wave period [s]
S.d50=0.0002;                                                              % KAMP & MILH & VR14 : median grain diameter [m]
S.porosity=0.4;                                                            % KAMP & MILH & VR14 : S.porosity (typically 0.4) [-]
S.tanbeta=0.05;                                                            % KAMP & MILH & VR14 : mean bed slope [ratio 1/slope]
S.reftime='2020-01-01';                                                    % Reference time (i.e. 'yyyy-mm-dd') <- leave empty to use t=0
S.endofsimulation='2022-01-01';                                            % time step [year] (use 1/365/4 for time-series) (not used in combination with adaptive time step)
S.ds0=100;                                                                 % initial space step [m]
S.Courant=1;
S.tc=0;                                                                    % switch for automatic time step
S.dt=2/365;                                                                 % time step [year]
S.struct=0;                                                                % switch for using hard structures
S.twopoints=2;                                                             % upwind treatment involving two points (1) or 1 point (0)
S.LDBstructures='';                                                        % LDB with hard structures ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
S.growth=0;                                                                % calibration of nourishment growth rate
S.nourish=0;                                                               % switch (0/1) for nourishments
S.LDBnourish='';                                                           % LDB with nourishment locations ([Nx2] ASCII FILE WITHOUT HEADER) (i.e. polygon around relevant grid cells) <- leave empty to use interactive mode!
S.nourrate=100;
S.nourstart=0;
S.spit_width=100;                                                          % width of tip of spit
S.xlimits=[-1000,1500];                                                    % X limits of plot [m] [1x2] <- leave empty to automatically do this                                                          % Y limits of plot [m] [1x2] <- leave empty to automatically do this
S.ylimits=[0,10000];
S.XYwave =[];                                                              % X,Y location and scale of wave arrow [m] [1x2] (automatically determined on the basis of xy-limits and wave angle
S.XYoffset=[0,0];                                                          % shift in X,Y locaton for plotting <- leave empty to automatically shift grid <- use [0,0] for no shift
S.pauselength=[];                                                          % pause between subsequent timesteps (e.g. 0.0001) <- leave empty to not pause plot
S.fignryear=4;
S.storageinterval=365/4;                                                   % Time interval of storage of output file ('output.mat'; [day])
S.outputdir='Output\spit_base\';                                           % output directory for plots and animation
S.debug=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN SHORELINES MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
S=ShorelineS(S);
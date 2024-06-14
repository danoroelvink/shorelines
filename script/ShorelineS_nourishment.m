%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..\functions\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;
S.LDBcoastline='Data\sandmotor11_2011-08-03open.xy';                       % LDB with initial coastline shape [Nx2] <- leave empty to use interactive mode!
S.Hso=1;                                                                   % wave height [m]
S.phiw0=310;                                                               % deep water wave angle [°N]
S.phif=[];                                                                 % Orientation of the foreshore [°N]
S.spread=10;                                                               % wave spreading [°] (wave_dir from range:  phiw0 +/- 0.5*spread)
S.randomseed=100;                                                          % make sure the same random series is generated
S.d=10;                                                                    % active profile height [m]
S.trform='CERC';                                                           % switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
S.b=.75e6;                                                                 % CERC : coeff in simple cerc formula
S.reftime='2011-08-01';                                                    % Reference time <- leave empty to use t=0
S.endofsimulation='2012-08-01';                                            % End time of simulation
S.tc=0.9;
%S.dt=3/24/365;
S.smoothfac=0;
S.ds0=100;                                                                 % initial space step [m]
S.twopoints=1;                                                             % upwind treatment involving two points (1) or 1 point (0)
S.spit_width=50;                                                           % width of tip of spit
S.xlimits=[68700,74900];                                                   % X limits of plot [m] <- leave empty to automatically do this
S.ylimits=[450000,455200];                                                 % Y limits of plot [m] <- leave empty to automatically do this
S.outputdir='output\nourishment\';                                         % output directory for plots and animation
S.storageinterval=365/4;                                                   % Time interval of storage of output file ('output.mat'; [day])
S.nourish=1;                                                               % switch (0/1) for nourishments
S.growth=1;                                                                % calibration of nourishment growth rate
S.LDBnourish='Data\zm_nour.nor';                                           % LDB with nourishment locations ([Nx2] ASCII FILE WITHOUT HEADER) (i.e. polygon around relevant grid cells) <- leave empty to use interactive mode!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN SHORELINES MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S,O]=ShorelineS(S);
save('output\Results.mat','S','O');

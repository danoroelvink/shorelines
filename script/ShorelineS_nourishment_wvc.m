%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..\functions\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;
S.LDBcoastline='Data\sandmotor11_2011-08-03open.xy';                       % LDB with initial coastline shape [Nx2] <- leave empty to use interactive mode!
S.WVCfile='Data\sandmotor.wvc';                                            % 
S.phif=[71385 450435 312; ...                                              % Orientation of the foreshore [°N]
        71920 451495 300; ...
        72360 452445 312; ...
        73008 453323 330; ... 
        73909 453362 312]; 
S.randomseed=100;                                                          % make sure the same random series is generated
S.d=10;                                                                    % active profile height [m]
S.ddeep=10;
S.dnearshore=6;
S.trform='CERC3';                                                          % switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
S.b=.75e6;                                                                 % CERC : coeff in simple cerc formula
S.qscal=0.5;                                                               % calibration coefficient for the transport
S.reftime='2011-08-01';                                                    % Reference time <- leave empty to use t=0
S.endofsimulation='2012-08-01';                                            % End time of simulation
S.tc=0;
S.dt=3/24/365;
S.smoothfac=0;
S.smoothrefrac=0.4;
S.ds0=100;                                                                 % initial space step [m]
S.twopoints=1;                                                             % upwind treatment involving two points (1) or 1 point (0)
S.spit_width=50;                                                           % width of tip of spit
S.xlimits=[68700,74900];                                                   % X limits of plot [m] <- leave empty to automatically do this
S.ylimits=[450000,455200];                                                 % Y limits of plot [m] <- leave empty to automatically do this
S.relaxationlength=2*S.ds0;                                                % length over which transport decelerates in meters, which adds inertia to the longshore current. It scales linearly with the wave height below 1m waves.
S.outputdir='output\sandmotor_wvc\';                                       % output directory for plots and animation
S.storageinterval=365/4;                                                   % Time interval of storage of output file ('output.mat'; [day])
S.plotinterval=1;
S.suppress_highangle=0;
S.plotUPW=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN SHORELINES MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S,O]=ShorelineS(S);
save('output\Results.mat','S','O');

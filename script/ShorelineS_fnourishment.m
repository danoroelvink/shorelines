%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..\functions\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;
S.LDBcoastline='data\fnour.txt';                                              % LDB with initial coastline shape [Nx2] <- leave empty to use interactive mode!
S.d=10;                                                                    % active profile height [m]
S.struct = 1;
S.x_hard = [5000 5000 5010 5010 nan 11000 11000 11010 11010 ];
S.y_hard = [-100 1000 1000 -100 nan -100 1000 1000 -100];

%% ----------------transport formula---------------------------------------
S.trform='CERC';                                                           % switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
S.qscal=1;                                                                 % Calibration factor of the transport (works for all transport formulas)

%% ----------------boundary conditions-------------------------------------
S.Hso = 2; 
S.phiw0 = 0; 
S.spread = 0; 
S.ddeep=10;                                                                % Waterdepth the location of wave climate, corresponding with S.Hso [m]
S.dnearshore=8;                                                            % Waterdepth at the 'dynamic boundary', corresponding with S.phif [m]

%% -------------------------- boundary condition -------------------------- for (non cyclic) sections (ex.straight shoreline) ,If more than one (non cyclic)sections should adjusted manually
S.boundary_condition_start='fixed';                                        % boundary condition 'CTAN', 'FIXD', 'FIXD2', 'FUNC'
S.boundary_condition_end='fixed';

%% -----------------time settings------------------------------------------
S.tc=0;                                                                    % switch for using adaptive time step: 0 (no adaptive time step), ~=0 (adaptive time step), 0.9 (factor on time step: close to 1 means larger time step)
S.dt=1/365;                                                                % time step (in days)
S.ds0=100;                                                                 % initial space step [m]
S.reftime='2000-01-01';                                                    % Reference time <- leave empty to use t=0
S.endofsimulation='2000-07-01';                                            % End time of simulation
S.twopoints=1;                                                             % upwind treatment involving two points (1) or 1 point (0)
S.smoothfac=0;

%% ------------------nourishment-------------------------------------------
S.nourish=0;                                                               % switch (0/1) for nourishments
S.LDBnourish = ''; 

%% ------------------shoreface nourishment---------------------------------
S.fnourish=1; 
S.fnorfile='data\fnour.fnor';
S.K=[0.7 0.7 0.7];

%% ------------output------------------------------------------------------
S.xlimits=[0 15000];                                                       % X limits of plot [m] <- leave empty to automatically do this
S.ylimits=[-2000 5000];                                                    % Y limits of plot [m] <- leave empty to automatically do this
S.outputdir='output\fnourishment\';                                        % output directory for plots and animation
S.storageinterval=30;                                                      % every 30 days

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN SHORELINES MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
[S,O]=ShorelineS(S);

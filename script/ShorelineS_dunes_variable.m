addpath(genpath('..\functions\'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;
S.LDBcoastline='Data\test_coastline_hor.xy';                               % LDB with initial coastline shape [Nx2] <- leave empty to use interactive mode!
S.WVCfile={'Data\wave_climate_015.wvc',-100,-250;'Data\wave_climate_015.wvc',-100,250};       % wave time-series <-leave empty to use wave parameters ('S.Hso', 'S.phiw0' and 'S.spread')
S.spread=0;                                                                % wave spreading [°] (wave_dir from range:  phiw0 +/- 0.5*spread)
S.phif=[];                                                                 % Orientation of the foreshore [°N]
S.d=8;                                                                     % active profile height [m]
S.trform='CERC';                                                           % switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
S.b=.75e6;                                                                 % CERC : coeff in simple cerc formula
S.reftime='2000-01-01';                                                    % Reference time <- leave empty to use t=0
S.endofsimulation='2022-05-30';                                            % End time of simulation
S.ddeep=5.58;                                                              % Waterdepth the location of wave climate, corresponding with S.Hso [m]
S.tc=0;                                                                    % adaptive time step
S.dt=1/365/2;                                                              % time step of the coastal processes in fraction of a year
S.dtdune=1/365/24;                                                         % time step of the dune processes in fraction of a year
S.smoothfac=0;
S.ds0=100;                                                                 % initial space step [m]
S.twopoints=1;                                                             % upwind treatment involving two points (1) or 1 point (0)
S.growth=0;
S.nourish=0;                                                               % switch (0/1) for nourishments
S.LDBnourish='';                                                           % LDB with nourishment locations [Nx2] (i.e. polygon around relevant grid cells) <- leave empty to use interactive mode!
S.nourrate=100;
S.nourstart=0;
S.spit_width=50;                                                           % width of tip of spit
S.xlimits=[-120 120];                                                      % X limits of plot [m] <- leave empty to automatically do this
S.ylimits=[-300 300];                                                      % Y limits of plot [m] <- leave empty to automatically do this
S.outputdir='Output\dunes_variable\';                                      % output directory for plots and animation
S.d50=2.5e-4;                                                              % KAMP & MILH & VR14 : median grain diameter [m]
S.boundary_condition_start='Fixed';                                        % boundary condition: 'Angleconstant', 'Fixed' , 'Closed', 'Periodic' , see "transport_boundary_condition.m"
S.boundary_condition_end='Fixed';

%% Extra inputs for dune evolution simulation
S.dune = 1; 
S.z=10;                                                                    % elevation of measred wind data
S.Cs=3e-3;                                                                 % Impact coefficient waves based on thesis M. Ghonim || very sensitive
S.LDBdune='Data\test_dunefoot_hor.dun';                                    % LDB with initial dune foot shape [Nx2] <- leave empty to use interactive mode!
S.WNDfile={'Data\wind_climate_015.wnd',-100,-250; 'Data\wind_climate_015.wnd',-100,250};                                     %'Windclimate\wind_climate.txt';       % wind time-series <-leave empty to use wave parameters ('S.uz', 'S.phiwnd0' )
S.WATfile={'Data\tide_015.wat',-100,-250; 'Data\tide_015.wat',-100,250};   % Water levels file relative to MSL
S.WVDfile={'Data\wave_climate_015.wvd',-100,-250; 'Data\wave_climate_015.wvd',-100,250};                                    % wave time-series <-leave empty to use wave parameters ('S.Hso', 'S.phiw0' and 'S.spread')
S.xyprofiles=[0,-200;0,0;0,200];
S.fignryear=4;
S.plotinterval=91;
S.storageinterval=91;                                                      % storage of output in days

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN SHORELINES MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
[S,O]=ShorelineS(S);

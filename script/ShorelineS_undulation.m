addpath(genpath('..\functions\'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;
S.Hso=1;                                                                   % wave height [m]
S.phiw0=300;                                                               % deep water wave angle [°N]
S.spread=10;                                                                % wave spreading [°] (wave_dir from range:  phiw0 +/- 0.5*spread)
S.tper=8;
S.reftime='2000-01-01';                                                    % Reference time <- leave empty to use t=0
S.endofsimulation='2002-01-02';                                            % time step [year] (use 1/365/4 for time-series) (not used in combination with adaptive time step)
S.LDBcoastline='Data\Undulation-0yrs.xy';                                  % LDB with initial coastline shape [Nx2] <- leave empty to use interactive mode!
S.boundary_condition_start='Fixed';                                     % boundary condition 'PRDC' for periodeic boundary condition or 'FIXD' for fixed boundary condition
S.boundary_condition_end='Fixed';
S.trform='CERC3'; 
S.ddeep=8;
S.dnearshore=4;
S.b=.5e6;                                                                  % CERC : coeff in simple cerc formula
S.d=1;                                                                     % active profile height [m]
S.ds0=250;
S.twopoints=1;                                                             % switch for 'twopoints appraoch'
S.spit_width=250;                                                          % width of tip of spit
S.xlimits=[25000,75000];                                                   % X limits of plot [m] <- leave empty to automatically do this
S.ylimits=[-5000,5000];                                                    % Y limits of plot [m] <- leave empty to automatically do this
S.XYwave = [];                                                             % X,Y location and scale of wave arrow [m] (automatically determined on the basis of xy-limits and wave angle
S.XYoffset=[0,0];
S.outputdir='Output\undulation\';
S.LDBplot ='';
S.plotinterval=1;
S.fignryear=1;         
S.plotUPW=1;
S.smoothrefrac=0;      
S.relaxationlength=500;                                                    % length over which transport decelerates in meters, which adds inertia to the longshore current. It scales linearly with the wave height below 1m waves.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN SHORELINES MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
S=ShorelineS(S);

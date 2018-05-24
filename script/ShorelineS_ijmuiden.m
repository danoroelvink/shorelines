
clear all;close all
addpath(genpath('ShorelineS\'))
addpath(genpath('..\functions\'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S.LDBcoastline='ijmuiden\xy1967_StN.txt';          
                                % LDB with initial coastline shape ([Nx2] ASCII FILE WITHOUT HEADER)
                                % <- leave empty to use interactive mode!
S.struct=1;                     % switch for using hard structures
S.LDBstructures='ijmuiden\ijmuiden_port.txt';       
                                % LDB with hard structures ([Nx2] ASCII FILE WITHOUT HEADER) 
                                % <- leave empty to use interactive mode!
S.LDBplot = {'ijmuiden\xy1967.txt','1967','k--';...
             'ijmuiden\xy2007.txt','2007','r--';...
             };                 % Filenames, dates and line codes of o=bserved coastlines
S.xlimits=[90000 110000];       % x-limits of plot area
S.ylimits=[490000 505000];      % y-limits of plot area
S.reftime='1967-07-01';         % Reference time <- leave empty to use t=0
S.endofsimulation='2007-07-01'; % End time of simulation
S.outputdir='Output_IJmuiden100';  % output directory for plots and animation
S.Waveclimfile='Waveclimate\wcon.txt';              
                                % Wave climate file <-leave empty to use wave parameters 
                                % ('S.Hso', 'S.phiw0' and 'S.spread')
S.phiw0=330;                    % deep water wave angle [°N]
S.spread=90;                    % wave spreading [°] (wave_dir from range:  S.phiw0 +/- 0.5*S.spread)
S.WVCfile='';                   % wave time-series <-leave empty to use wave parameters ('S.Hso', 'S.phiw0' and 'S.spread')
S.d=10;                         % Active profile height [m]
S.ddeep=30;                     % Offshore water depth for refraction
S.dnearshore=6;                 % Nearshore water depth for refraction
S.phif=[];                      % Orientation of the foreshore [°N] <- only relevant for 'KAMP', 'MILH' or 'VR14'
S.trform='CERC';               % Switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
S.b=0.5e6;                       % CERC : coeff in simple cerc formula
S.tc=1;                       % factor to modify automatic timestep
S.ds0=200;                      % initial space step [m]
S.nt=5000;                      % number of timesteps
S.growth=0;                     % calibration of nourishment growth rate
S.nourish=0;                    % switch (0/1) for nourishments
S.LDBnourish='';                % LDB with nourishment locations ([Nx2] ASCII FILE WITHOUT HEADER) (i.e. polygon around relevant grid cells) <- leave empty to use interactive mode!
S.nourrate=100;                 % nourishment rate in polygon (m3/m/y)
S.nourstart=0;                  % start time of nourishment (y)
S.spit_width=50;                % width of tip of spit
S.twopoints=0;                  % upwind treatment involving two points (1) or 1 point (0)
S.plotinterval=1;
S.fignryear=1;
S.smoothfac=0.1;

S.phiw0=S.phiw0*pi/180;
S.phif=S.phif*pi/180;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN SHORELINES MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S]=ShorelineSm4(S);
save('Results','S')
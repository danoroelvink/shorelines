clear all;close all
addpath(genpath('ShorelineS\'))
addpath(genpath('..\functions\'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S.LDBcoastline='Data\LDB_161\1984-01-09-1984-07-07_IM_L4_UTM_SLini_South.txt';          
                                % LDB with initial coastline shape ([Nx2] ASCII FILE WITHOUT HEADER)
                                % <- leave empty to use interactive mode!
S.start=0;
S.struct=1;                     % switch for using hard structures
S.channel=1;
S.LDBchannel='Data\Senegal.txt';
S.channel_width=[300];
S.channel_fac=0.2;
S.x_hard=[337200 337200 350000];
S.y_hard=[1754600 1755100 1755100];
% if ~S.start>0
%    load rivers.mat;
% end
% S.xr_mc=xr_mc;
% S.yr_mc=yr_mc;
                                % LDB with hard structures ([Nx2] ASCII FILE WITHOUT HEADER) 
                                % <- leave empty to use interactive mode!
S.wavetransm = [0,0];           % transmission coeficient of the structures
                                % array size equal to number of structure
                                % sections
S.LDBplot = {'Data\LDB_161\1984-01-09-1984-07-07_IM_L4_UTM.txt','1984','k:';...
             'Data\LDB_161\1984-10-23-1985-04-21_IM_L4_UTM.txt','1985','b:';...
             'Data\LDB_161\1988-08-31-1989-02-27_IM_L4_UTM.txt','1988','g:';...
             'Data\LDB_161\2003-02-06-2003-08-05_IM_L4_UTM.txt','2003','r:';...
             };                 % Filenames, dates and line codes of o=bserved coastlines
S.xlimits=[ 330000  345000];    % x-limits of plot area
S.ylimits=[1745000 1755000];    % y-limits of plot area
S.reftime='1984-04-09';         % Reference time <- leave empty to use t=0
S.endofsimulation='2003-04-01'; % End time of simulation
S.outputdir='StLouisCal310';       % output directory for plots and animation
S.Waveclimfile='';              
                                % Wave climate file <-leave empty to use wave parameters 
                                % ('S.Hso', 'S.phiw0' and 'S.spread')
S.Hso_in=1.4;                     % deep water wave height(m)
S.Hso=1.4;                        % deep water wave height(m)
S.phiw0=310;                    % deep water wave angle [°N]
S.spread=60;                    % wave spreading [°] (wave_dir from range:  S.phiw0 +/- 0.5*S.spread)
S.WVCfile='';                   % wave time-series <-leave empty to use wave parameters ('S.Hso', 'S.phiw0' and 'S.spread')
S.d=5;                         % Active profile height [m]
S.ddeep=30;                     % Offshore water depth for refraction
S.dnearshore=6;                 % Nearshore water depth for refraction
S.phif=[];                      % Orientation of the foreshore [°N] <- only relevant for 'KAMP', 'MILH' or 'VR14'
S.trform='CERC';               % Switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
S.b=.5e6;                        % CERC : coeff in simple cerc formula
S.tc=.9;                       % factor to modify automatic timestep
S.ds0=100;                      % initial space step [m]
S.nt=5000;                      % number of timesteps
S.growth=0;                     % calibration of nourishment growth rate
S.nourish=0;                    % switch (0/1) for nourishments
S.LDBnourish='';                % LDB with nourishment locations ([Nx2] ASCII FILE WITHOUT HEADER) (i.e. polygon around relevant grid cells) <- leave empty to use interactive mode!
S.nourrate=100;                 % nourishment rate in polygon (m3/m/y)
S.nourstart=0;                  % start time of nourishment (y)
S.spit_width=1;                 % width of tip of spit
S.twopoints=2;                  % upwind treatment involving two points (1) or 1 point (0)
S.plotinterval=10;
S.fignryear=12;
S.smoothfac=0.05;
S.usefill=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN SHORELINES MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S]=ShorelineS(S);
save('Results','S')
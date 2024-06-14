%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('..\functions\'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;

%% Parameters that are varied over the 16 tests
%% Fixed parameters
S.yesplot=0;
S.x_mc=[0 2500];
S.y_mc=[0 0];
S.x_hard=[1100 1400 ];
S.y_hard=[250 250 ];
S.reftime='2022-01-01';
S.endofsimulation='2022-02-01';
S.Hso=1;                                                                   % wave height [m]
S.phiw0=40;                                                                % deep water wave angle [�N]
S.phif=360;                                                                % orientation of the deeper shoreface [degrees North]
S.ddeep=8;                                                                 % offshore depth [m]
S.dnearshore=4;                                                            % nearshore depth [m]
S.d=2;                                                                     % active profile height [m] <- this variable is low only for demonstration purposes!
S.b=.5e6;                                                                  % CERC : coeff in simple cerc formula%
S.tc=0;                                                                    % adaptive time step CFL
S.dt=3/24/365;
S.boundary_condition_start='Neumann';
S.boundary_condition_end='Neumann';
S.ds0=50;                                                                  % initial space step [m]
S.struct=1;                                                                % switch for using hard structures
S.twopoints=1;
S.xlimits=[0 2500];
S.ylimits=[-250 500];
S.plotinterval=1;
S.storageinterval=365/8;
S.fignryear=36;
S.smoothfac=0.0;
S.ld=1000;
S.video=0;
S.spit_width=50;
S.randomseed=100;                                                          % make sure the same random series is generated
S.usefill=1;
S.spread=0;
S.Aw=5;
S.trform='CERC3';
S.diffraction=1;
S.outputdir=['output\offshorebreakwater\'];
S.storageinterval=10;                                                      % Time interval of storage of output file ('output.mat'; [day])
S.plotHS=6;
S.plotDIR=6;
S.plotQS=0;
S.plotUPW=0;
S.xyout=[0,250,2000,250];
S.relaxationlength=2*S.ds0;                                                % length over which transport decelerates in meters, which adds inertia to the longshore current. It scales linearly with the wave height below 1m waves.
S.maxangle=60;

%% RUN SHORELINES MODEL
[S,O]=ShorelineS(S);

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
S.x_hard=[900 1600 ];
S.y_hard=[300 300 ];
S.reftime='2022-01-01';
S.endofsimulation='2023-01-01';
S.Hso=1;                                                                   % wave height [m]
S.phiw0=360;                                                                 % deep water wave angle [�N]
S.spread=0;
S.phif=360;                                                                % orientation of the deeper shoreface [degrees North]
S.ddeep=20;                                                                % offshore depth [m]
S.dnearshore=6;                                                            % nearshore depth [m]
S.d=6;                                                                     % active profile height [m] <- this variable is low only for demonstration purposes!
S.b=.5e6;                                                                  % CERC : coeff in simple cerc formula%
S.tc=0;                                                                    % adaptive time step CFL
S.dt=6/24/365;
S.boundary_condition_start='Closed';
S.boundary_condition_end='Closed';
S.ds0=50;                                                                  % initial space step [m]
S.struct=1;                                                                % switch for using hard structures
S.twopoints=1;
S.xlimits=[0 2500];
S.ylimits=[-150 550];
S.plotinterval=1;
S.storageinterval=365/8;
S.fignryear=36;
S.smoothfac=0.0;
S.ld=1000;
S.video=0;
S.spit_width=50;
S.usefill=1;
S.dirspr=12;
S.randomseed=100;                                                          % make sure the same random series is generated
S.Aw=5;
S.trform='CERC3';
S.diffraction=1;
S.outputdir=['output\offshorebreakwater_wavespread\'];
S.storageinterval=30;                                                      % Time interval of storage of output file ('output.mat'; [day])
S.plotHS=0;
S.plotDIR=3;
S.plotQS=0;
S.plotUPW=0;

%% RUN SHORELINES MODEL
[S,O]=ShorelineS(S);

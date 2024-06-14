%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('..\functions\'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;
S.reftime='2020-01-01';                                                    % Reference time (i.e. 'yyyy-mm-dd') <- leave empty to use t=0
S.endofsimulation='2021-01-01';
S.LDBcoastline='Data\circle.txt';
S.Waveclimfile='Data\wcon.txt';
S.randomseed=100;                                                          % make sure the same random series is generated
S.d=3;                                                                     % active profile height [m] <- the value used here is very small for the purpose of quickly showing a demonstration
S.b=1e6;                                                                   % CERC : coeff in simple cerc formula
S.tc=0;                                                                    % time step fator
S.dt=1/365;
S.ds0=200;                                                                 % initial space step [m]
S.struct=0;                                                                % switch for using hard structures
S.spit_width=100;                                                          % minimum width for overwash
S.xlimits=[0 10000];
S.ylimits=[0 10000];
S.twopoints=1;
S.outputdir='Output\circle\';                                              % output directory for plots and animation
S.plotinterval=10;
S.fignryear=8;
S.storageinterval=30;                                                      % Time interval of storage of output file ('output.mat'; [day])
S.smoothfac=0.01;
S.plotvisible=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN SHORELINES MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S,O]=ShorelineS(S);

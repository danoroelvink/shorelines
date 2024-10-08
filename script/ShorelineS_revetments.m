addpath(genpath('..\functions\'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;
S.reftime='2017-01-01';
S.endofsimulation='2017-05-01';
S.phif=270;
S.d=6;                                                                     % active profile height [m]
S.ddeep=10;
S.dnearshore=10;
S.b=.5e6;                                                                  % CERC : coeff in simple cerc formula%
S.tc=0;                                                                    % adaptive time step
S.dt=1/365/2;
S.boundary_condition_start='Neumann';
S.boundary_condition_end='Neumann';
S.griddingmethod=1;
S.ds0=50;                                                                  % initial space step [m]
S.spit_width=100;
S.OWscale=0.01;
S.trform='CERC3';
S.revet=1;   % switch for using hard structures
S.crit_width=5;
S.diffraction=1;
S.twopoints=1;
S.xlimits=[0 1000];
S.ylimits=[700 2300];
S.plotinterval=10;
S.storageinterval=20;
S.fignryear=12;
S.smoothfac=0;
S.ld=1000;
S.video=0;
S.usefill=1;
S.plotDIR=4;
S.plotUPW=1;

%% Test 1
test=1;
S.x_mc=[500 500 1000];
S.y_mc=[0 1500 3000];
S.x_revet=[550 550 650];
S.y_revet=[1200 1500 1800];
S.Hso=0.7;                                                                 % wave height [m]
S.phiw0=310;                                                               % deep water wave angle [�N]
S.spread=40;                                                               % wave spreading [�] (wave_dir from range:  phiw0 +/- 0.5*spread)
S.randomseed=100;
S.outputdir=['output\revetments',num2str(test,'%3.3i')];
[S]=ShorelineS(S); % RUN SHORELINES MODEL

%% Test 2
test=2;
S.x_mc=[1000 500 1000];
S.y_mc=[0 1500 3000];
S.x_revet=[650 550 650];
S.y_revet=[1200 1500 1800];
S.Hso=0.7;                                                                 % wave height [m]
S.phiw0=270;                                                               % deep water wave angle [�N]
S.spread=40;                                                               % wave spreading [�] (wave_dir from range:  phiw0 +/- 0.5*spread)
S.randomseed=100;
S.spit_width=1;
S.outputdir=['output\revetments',num2str(test,'%3.3i')];
[S]=ShorelineS(S); % RUN SHORELINES MODEL

%% Test 3
test=3;
S.x_mc=[500 1000 500];
S.y_mc=[0 1500 3000];
S.x_revet=[850 950 850];
S.y_revet=[1200 1500 1800];
S.Hso=0.7;                                                                 % wave height [m]
S.phiw0=270;                                                               % deep water wave angle [�N]
S.spread=0;                                                                % wave spreading [�] (wave_dir from range:  phiw0 +/- 0.5*spread)
S.outputdir=['output\revetments',num2str(test,'%3.3i')];
[S]=ShorelineS(S); % RUN SHORELINES MODEL

%% Test 4
test=4;
S.x_mc=[1000 500 1000];
S.y_mc=[0 1500 3000];
S.x_revet=[650 550 650]-50;
S.y_revet=[1200 1500 1800];
S.Hso=0.7;                                                                 % wave height [m]
S.phiw0=270;                                                               % deep water wave angle [�N]
S.spread=0;                                                                % wave spreading [�] (wave_dir from range:  phiw0 +/- 0.5*spread)
S.spit_width=1;
S.outputdir=['output\revetments',num2str(test,'%3.3i')];
[S]=ShorelineS(S); % RUN SHORELINES MODEL

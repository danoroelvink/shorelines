%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('..\functions\'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;

%% Parameters that are varied over the 16 tests
%% Fixed parameters
S.yesplot=0;
S.x_mc=[500 2500];
S.y_mc=[0 0];
S.x_hard=[490 490  700 700 510 510, nan , 1250 1750 , nan , 2490 2490 2300 2300 2510 2510 ];
S.y_hard=[-150 350 350 330 330 -150, nan, 350 350 , nan, -150 330 330 350 350 -150 ];
S.reftime='2022-01-01';
S.endofsimulation='2022-04-01';
S.Hso=1;                                                                   % wave height [m]
S.phiw0=0;                                                                 % deep water wave angle [�N]
S.phif=360;                                                                % orientation of the deeper shoreface [degrees North]
S.ddeep=20;                                                                % active profile height [m]
S.dnearshore=6;                                                            % active profile height [m]
S.d=6;                                                                     % active profile height [m]
S.b=.5e6;                                                                  % CERC : coeff in simple cerc formula%
S.tc=0;                                                                    % adaptive time step CFL
S.dt=6/24/365;
S.boundary_condition_start='Closed';
S.boundary_condition_end='Closed';
S.griddingmethod=2;
S.ds0=50;                                                                  % initial space step [m]
S.struct=1;                                                                % switch for using hard structures
S.twopoints=1;
S.xlimits=[0 3000];
S.ylimits=[-250 500];
S.plotinterval=1;
S.storageinterval=30;                                                      % every 30 days
S.fignryear=36;
S.smoothfac=0.0;
S.ld=1000;
S.video=0;
S.spit_width=50;
S.usefill=1;
S.spread=50;
S.Aw=5;
S.trform='CERC3';
S.diffraction=1;
S.outputdir='output\groynes2obw_tshape_bnd\';
S.plotHS=6;
S.plotDIR=6;
S.plotQS=0;
S.plotUPW=0;

%% RUN SHORELINES MODEL
[S,O]=ShorelineS(S);

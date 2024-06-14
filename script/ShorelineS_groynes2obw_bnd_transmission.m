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
S.x_hard=[490 490 510 510, nan , 1250 1750 , nan , 2490 2490 2510 2510 ];
S.y_hard=[-150 350 350 -150, nan, 350 350 , nan, -150 350 350 -150 ];
S.transmission=1;                                                          % this switch controls whether breakwaters allow for waves to propagate over them
S.transmform='angr';                                                       % Wave transmission approach 'angr_perm', 'angr_imperm', 'gent_perm', 'gent_imperm', 'seabrhall'          
S.transmdir=1;                                                             % switch to use original wave direction after transmission over breakwater [S.BW_origdir=1], or breakwater orientation [S.BW_origdir=0];
S.transmbwdepth=[8];                                                       % if S.transmission=1 - this defines at which depth breakwater is positioned [m] as input for for wave_transmission.m function Positive downward. 
S.transmcrestheight=-0.75;                                                 % if S.transmission=1 - this defines breakwater height [m] w.r.t. MSL, as input for for wave_transmission.m function. Positive upward, above MSL.
S.transmslope=1/3;                                                         % if S.transmission=1 - this defines breakwater slope [-] as input for for wave_transmission.m function 
S.transmcrestwidth=6;                                                      % if S.transmission=1 - this defines breakwater width [m] as input for for wave_transmission.m function
S.WATfile='Data\tide_015.wat';                                             % Water levels file relative to MSL

S.reftime='2022-01-01';
S.endofsimulation='2022-04-01';
S.Hso=1;                                                                   % wave height [m]
S.phiw0=20;                                                                % deep water wave angle [�N]
S.phif=360;                                                                % orientation of the deeper shoreface [degrees North]
S.ddeep=20;                                                                % active profile height [m]
S.dnearshore=6;                                                            % active profile height [m]
S.d=6;                                                                     % active profile height [m]
S.b=.5e6;                                                                  % CERC : coeff in simple cerc formula%
S.tc=0;                                                                    % adaptive time step CFL
S.dt=6/24/365;
S.boundary_condition_start='Closed';
S.boundary_condition_end='Closed';
S.griddingmethod=1;
S.ds0=50;                                                                  % initial space step [m]
S.struct=1;                                                                % switch for using hard structures
S.twopoints=1;
S.xlimits=[0 3000];
S.ylimits=[-250 500];
S.plotinterval=20;
S.storageinterval=30;                                                      % every 30 days
S.fignryear=36;
S.smoothfac=0.0;
S.ld=1000;
S.video=0;
S.spit_width=50;
S.usefill=1;
S.spread=0;
S.Aw=5;
S.trform='CERC3';
S.diffraction=1;
S.outputdir='output\groynes2obw_bnd_transmission\';
S.plotHS=6;
S.plotDIR=6;
S.plotQS=0;
S.plotUPW=0;

%% RUN SHORELINES MODEL
[S,O]=ShorelineS(S);

addpath(genpath('..\functions\'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;
S.reftime='2020-01-01';
S.endofsimulation='2030-01-01';
S.Hso=1;                                                                   % wave height [m]
S.phiw0=270;                                                               % deep water wave angle [°N]
S.phif=270;
S.spread=0;                                                               % wave spreading [°] (wave_dir from range:  phiw0 +/- 0.5*spread)
S.d=6;                                                                     % active profile height [m]
S.b=.5e6;                                                                  % CERC : coeff in simple cerc formula%
S.tc=0.9;                                                                    % adaptive time step
%S.dt=1/365;                                                                % time step [yr]
S.ds0=200;                                                                 % initial space step [m]
S.struct=1;                                                                % switch for using hard structures
S.LDBstructures=[];
S.x_mc=[4000 4000];
S.y_mc=[0 5000];
S.x_hard=3400-500*cosd([0:30:360]);
S.y_hard=2500+500*sind([0:30:360]);
S.diffraction=1;
S.twopoints=1;
S.spit_width=50;                                                           % width of tip of spit
S.xlimits=[0 8000];
S.ylimits=[0 5000];
S.outputdir='Output\interactive\';
S.plotinterval=1;
S.storageinterval=1e6;
S.fignryear=1;
S.smoothfac=0.0;
S.ld=1000;
S.video=0;
S.usefill = 0;
S.debug=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN SHORELINES MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S]=ShorelineS(S);

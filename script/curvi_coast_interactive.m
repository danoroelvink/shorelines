addpath(genpath('ShorelineS\'))
addpath(genpath('..\functions\'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;
S.reftime='2017-04-09';
S.endofsimulation='2067-04-09';
S.Hso=1;                                                                   % wave height [m]
S.phiw0=270*pi/180;                                                        % deep water wave angle [°N]
S.spread=360;                                                               % wave spreading [°] (wave_dir from range:  phiw0 +/- 0.5*spread)
S.d=10;                                                                    % active profile height [m]
S.b=.5e6;                                                                   % CERC : coeff in simple cerc formula
S.dt=1/50;                                                                  % time step [year]
S.ds0=50;                                                                 % initial space step [m]
S.nt=5000;                                                                 % number of timesteps
S.struct=1;                                                                % switch for using hard structures
S.twopoints=1;
S.spit_width=50; % width of tip of spit
S.xlimits=[0 3000];
S.ylimits=[0 3000];
ntt=1;
S.seaslope=0.01;
S.landslope=0.05;
S.seamin=-3;
S.landmax=3;
S.tide_interaction=false;
S.dx=200;
S.dy=200;
S.Lx=10000;
S.Ly=40000;
S.zdeep=-10;
S.zshallow=-2;
S.slope=.002;
S.outputdir='Output_interactive';
S.plotinterval=1;
S.fignryear=1;
S.smoothfac=0.1;
S.ld=1000;
S.video=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN SHORELINES MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xg,yg,zg]=inigrid( S.dx,S.dy,S.Lx,S.Ly,S.zdeep,S.slope,S.zshallow );

for itt=1:ntt
    [S]=ShorelineSm4(S);
    if S.tide_interaction
        [zg]=update_bathy(xg,yg,zg,S);
        !run
        [S,zg]=update_shoreline(xg,yg,zg,S);
    end
end
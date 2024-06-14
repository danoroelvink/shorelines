%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('..\functions\'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;

%% Parameters that are varied over the 16 tests
%% Fixed parameters
S.yesplot=0;
S.x_mc=[0 3000];
S.y_mc=[0 0];
if 1
    S.x_hard=[1090 1090 1110 1110, nan , 2090 2090 2110 2110 ];
    S.y_hard=[-150 350 350 -150 , nan, -150 350 350 -150 ];
else
    S.x_hard=[1090 1090 1110 1110, nan , 2090, 2090, 1900, 1000, 1000, 1900, 2110, 2110 ];
    S.y_hard=[-150 350 350 -150,nan, -150,  350,  450,  450,  470,  470,  360, -150 ];
end

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
S.boundary_condition_start='Neumann';
S.boundary_condition_end='Neumann';
if 1 % regular grid
    S.griddingmethod=2;
    S.ds0=50;                                                              % initial space step [m]
else % use spatially varying grid (not used now)
    S.griddingmethod=1;
    S.ds0=[0,0,300; ...
           1300,0,50; ...
           2500,0,300];
end
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
S.spread=360;
S.randomseed=100;
S.Aw=5;
S.trform='CERC3';
S.diffraction=1;
S.outputdir='output\groynes2\';
S.plotHS=6;
S.plotDIR=6;
S.plotQS=0;
S.plotUPW=0;

%% RUN SHORELINES MODEL
[S,O]=ShorelineS(S);

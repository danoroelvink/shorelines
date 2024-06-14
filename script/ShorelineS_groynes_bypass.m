%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('..\functions\'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;

%% Parameters that are varied over the 16 tests
%% Fixed parameters
for ii=1:3
    S.yesplot=0;
    S.x_mc=[800 1700];
    S.y_mc=[0 0];
    S.x_hard=[1290 1290 1310 1310 ];
    S.y_hard=[-150 100 100 -150 ];
    S.reftime='2022-01-01';
    S.endofsimulation='2022-04-01';
    S.Hso=1;                                                                   % wave height [m]
    S.phiw0=40;                                                                % deep water wave angle [�N]
    S.phiw0=320;                                                               % deep water wave angle [�N]
    S.phif=360;                                                                % orientation of the deeper shoreface [degrees North]
    S.ddeep=20;                                                                % active profile height [m]
    S.dnearshore=6;                                                            % active profile height [m]
    S.d=6;                                                                     % active profile height [m]
    S.b=.5e6;                                                                  % CERC : coeff in simple cerc formula%
    S.tc=0;                                                                    % adaptive time step CFL
    S.dt=12/24/365;
    S.boundary_condition_start='Fixed';
    S.boundary_condition_end='Fixed';
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
    S.xlimits=[800 1800];
    S.ylimits=[-200 300];
    S.plotinterval=1;
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
    if ii==1
        S.bypasscontractionfactor=1;
    elseif ii==2
        S.bypasscontractionfactor=1.1;
    else
        S.bypasscontractionfactor=1.2;
    end
    S.outputdir=['output\groynes_bypass',num2str(S.bypasscontractionfactor,'%2.1f'),'\'];
    S.plotHS=6;
    S.plotDIR=6;
    S.plotQS=0;
    S.plotUPW=0;
    S.Aw=1.3; %1.3;
    S.Awfixedhs=5;
    %% RUN SHORELINES MODEL
    [S,O]=ShorelineS(S);
end

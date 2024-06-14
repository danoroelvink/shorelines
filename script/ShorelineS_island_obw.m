addpath(genpath('..\functions\'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ispw=1
    for ispr=2
        S=struct;
        S.x_hard=[280,300,300,280,280];
        S.y_hard=[-500,-500,500,500,-500];
        S.x_mc=[800,800];
        S.y_mc=[-4000,4000];
        S.reftime='2020-01-01';
        S.endofsimulation='2023-01-01';
        S.Hso=1;                                                                   % offshore significant wave height [m]
        S.phiw0=270;                                                               % deep water wave angle [°N]
        S.spread=(ispr-1)*90;                                                      % wave spreading [°] (wave_dir from range:  phiw0 +/- 0.5*spread)
        S.randomseed=100;                                                          % make sure the same random series is generated
        S.d=10;                                                                    % active profile height [m]
        S.b=1e6;                                                                   % CERC : coeff in simple cerc formula
        S.tc=0;                                                                    % switch for automatic time step
        S.dt=4/365;                                                                % time step [year]
        S.ds0=100;                                                                 % initial space step [m]
        S.struct=1;                                                                % switch for using hard structures
        S.twopoints=1;
        S.xlimits=[-1000 1500];
        S.ylimits=[-1500 1500];
        S.outputdir=['Output\island_obw',num2str(ispw),num2str(ispr),'\'];                                                  % output directory for plots and animation
        S.storageinterval=365/2;                                                   % Time interval of storage of output file ('output.mat'; [day])
        S.plotinterval=1;
        S.fignryear=4;
        S.smoothfac=0;
        S.spit_width=(ispw)*100;
        S.OWscale=0.02;
        S.debug=0;
        S.diffraction=1;
        S.maxangle=60;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% RUN SHORELINES MODEL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [S,O]=ShorelineS(S);
    end
end
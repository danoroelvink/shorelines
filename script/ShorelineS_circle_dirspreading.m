addpath(genpath('..\functions\'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ispw=1
    for ispr=1
        S=struct;
        for iang=1:73
            ang=-5*(iang-1);
            S.x_mc(iang)=0+500*cosd(ang);
            S.y_mc(iang)=0+500*sind(ang);
        end
        S.reftime='2020-01-01';
        S.endofsimulation='2021-01-01';
        S.Hso=1;                                                                   % offshore significant wave height [m]
        S.phiw0=270;                                                               % deep water wave angle [°N]
        S.spread=(ispr-1)*90;                                                      % wave spreading [°] (wave_dir from range:  S.phiw0 +/- 0.5*S.spread)
        S.randomseed=100;                                                          % make sure the same random series is generated
        S.d=10;                                                                    % active profile height [m]
        S.b=1e6;                                                                   % CERC : coeff in simple cerc formula
        S.tc=0;                                                                    % switch for automatic time step
        S.dt=0.5/365;                                                              % time step [year]
        S.ds0=50;                                                                 % initial space step [m]
        S.struct=0;                                                                % switch for using hard structures
        S.twopoints=1;
        S.xlimits=[-1000 1400];
        S.ylimits=[-1200 1200];
        S.outputdir=['Output\circle_dirspreading\'];
        S.storageinterval=30;                                                      % Time interval of storage of output file ('output.mat'; [day])
        S.plotinterval=10;
        S.fignryear=8;
        S.smoothfac=0;
        S.smoothrefrac=.5;
        S.griddingmethod=1;
        S.relaxationlength=80;
        S.spit_width=(ispw-1)*100+20;
        S.maxangle=60;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% RUN SHORELINES MODEL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [S]=ShorelineS(S);
    end
end

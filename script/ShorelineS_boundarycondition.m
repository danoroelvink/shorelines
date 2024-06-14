addpath(genpath('..\functions\'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=struct;
for mm=1:5
    if mm==1
        S.boundary_condition_start='Closed';
        S.boundary_condition_end='Angleconstant';
    elseif mm==2
        S.boundary_condition_start='Closed';
        S.boundary_condition_end={'Closed',-2e5};
    elseif mm==3
        S.boundary_condition_start='Angleconstant';
        S.boundary_condition_end='Fixed';
    elseif mm==4
        S.boundary_condition_start='Fixed';
        S.boundary_condition_end='Fixed';
    elseif mm==5
        S.boundary_condition_start='Closed';
        S.boundary_condition_end={'Angleconstant',20};
    end      
    S.x_mc=[0,2000];
    S.y_mc=[100,100];
    S.reftime='2020-01-01';
    S.endofsimulation='2021-01-01';
    S.Hso=1;                                                                   % wave height [m]
    S.phiw0=345;                                                               % deep water wave angle [°N]
    S.spread=0;
    S.d=2;                                                                     % active profile height [m] <- the value used here is very small for the purpose of quickly showing a demonstration
    S.b=1e6;                                                                   % CERC : coeff in simple cerc formula
    S.tc=0;                                                                    % switch for using adaptive time step
    S.ds0=100;                                                                 % initial space step [m]
    S.nt=5000;                                                                 % number of timesteps
    S.dt=1/365;                                                                % time step [year] -> use automatic timestep if S.dt==0 || S.tc==1
    S.struct=0;                                                                % switch for using hard structures
    S.twopoints=1;
    S.channel_width=0; 
    S.xlimits=[-50 2050];
    S.ylimits=[-100 500];
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
    if ischar(S.boundary_condition_end)
        S.outputdir=['Output\boundaryconditions_',num2str(mm,'%03.0f'),'_',lower(S.boundary_condition_start),'-',lower(S.boundary_condition_end)];
    else
        S.outputdir=['Output\boundaryconditions_',num2str(mm,'%03.0f'),'_',lower(S.boundary_condition_start),'-',lower(S.boundary_condition_end{1}),num2str(S.boundary_condition_end{2},'%+1.0f')];
    end
    S.storageinterval=30;                                                      % Time interval of storage of output file ('output.mat'; [day])
    S.plotinterval=60;
    S.fignryear=12;
    S.smoothfac=0.05;
    S.ld=1000;
    S.video=0;
    S.growth=0;
    S.spit_width=100;
    S.plotDIR=5;
    S.plotQS=5;
    S.xlimits=[0,2000];                                                   % X limits of plot [m] <- leave empty to automatically do this
    S.ylimits=[-300,500];                                                 % Y limits of plot [m] <- leave empty to automatically do this

    [S]=ShorelineS(S);
end

clear all
close all
clc

i=0
slos={'CERC','CERC2','KAMP','CERC3'};
blos={'CTAN','FIXD','FIXD2','PRDC','PRDC2'};
tc=[0.5,0.75,1,1.5];
phiwb=[0,-5,10]
phiwc=[0,5,15]
ypdo=[0:20:100]
xhrd=[50  75  100 125] 
yhrd=[ 150 200]
for iss=3
for izz=1
for iyy=3
for ioo=1
for ixh=4
for iyh=1
try
warning('Index exceeds matrix dimensions.');
 %webread('Index exceeds matrix dimensions.')

i=i+1;
addpath(genpath('ShorelineS\'))
addpath(genpath('..\functions\'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

S=struct;
S.Hso=1.2;% wave height [m]
S.WVCfile='';%wavedata_1967-2007.txt';
S.Waveclimfile='';
S.Wavecorr=0;%phiwc(izz); % Correction wave climate
S.RWSfiletype=true;
S.phiw0=phiwb(izz)*pi/180;                                                        % deep water wave angle [°N]
S.phif=S.phiw0;
S.spread=phiwc(izz);                                                               % wave spreading [°] (wave_dir from range:  phiw0 +/- 0.5*spread)
S.trform=slos{iss}%'CERC'; 
S.tper=7;                                                                  % KAMP & MILH : peak wave period [s]
S.d=10;                                                                    % active profile height [m]
% S.d=2.28*S.Hso-68.5*(S.Hso^2/9.81/S.tper^2);
S.LDBcoastline='s_shore800m.xy';                           % LDB with initial coastline shape [Nx2] <- leave empty to use interactive mode!
S.b=0.5e6;                                                                   % CERC : coeff in simple cerc formula
S.crit=.9;                                                                 % stability criterion (not active)
S.reftime='1970-05-01';                                                    % Reference time <- leave empty to use t=0
%S.dt=ddt(idt)/360;                                                                  % time step [year]
S.ds0=25;
S.tc=0.5;%tc(izz);                                                                          % number of timesteps
S.twopoints=0;                                                             % switch for 'twopoints appraoch'
S.struct=1;                                                                % switch for using hard structures
S.x_hard=[-xhrd(ixh) xhrd(ixh) ] ;
S.y_hard=[yhrd(iyh) yhrd(iyh) ];
S.LDBstructures='single_BW2.xy';                            % LDB with hard structures [Nx2] <- leave empty to use interactive mode!
S.growth=0;
S.nourish=0;                                                               % switch (0/1) for nourishments
S.LDBnourish='';                                                           % LDB with nourishment locations [Nx2] (i.e. polygon around relevant grid cells) <- leave empty to use interactive mode!
S.nourrate=100;
S.nourstart=0;
S.spit_width=20;                                                           % width of tip of spit
S.smoothfac=0.1;
S.xlimits=[-200,200];                                                  % X limits of plot [m] <- leave empty to automatically do this
S.ylimits=[-50,300];                                                % Y limits of plot [m] <- leave empty to automatically do this
S.XYwave = [];                                                             % X,Y location and scale of wave arrow [m] (automatically determined on the basis of xy-limits and wave angle
S.XYoffset=[0,0];
S.pauselength=[]; %0.0001;                                                 % pause between subsequent timesteps <- leave empty to not pause plot
S.outputdir=(strcat('Output\','_',S.trform));%,'_',num2str(k*100),'_ds0_',num2str(S.ds0),'_phi_',num2str(S.phiw0/pi*180)));
S.LDBplot ={};
S.plotinterval=3/360;
S.fignryear=196;
S.endofsimulation='1973-05-02';
S.BCfile='';
%% boundary condition                                                       for (non cyclic) sections (ex.straight shoreline) ,If more than one (non cyclic)sections should adjusted manually 
S.boundary_condition_start='FIXD2';%blos{izz};                                               % boundary condition 'PRDC' for periodeic boundary condition or 'FIXD' for fixed boundary condition
S.boundary_condition_end='FIXD2';%{izz};
%% Extract shorelines
S.SLplot={'1970-08-01','3 m','k--'
 '1970-11-01','6 m','b-'
 '1971-04-29','1 yr','g-'
 '1972-04-29','2 yr','K-'
 '1973-04-29','3 yr','r-'};
S.extract_x_y=0;                                                            %get file with shorelines coordinates x,y
S.print_fig=1;
%% video
S.video=1;
S.ypd=ypdo(ioo);
catch ME
%fprintf('WEBREAD without success: %s\n', ME.message);
    
continue

end

   
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN SHORELINES MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
S=ShorelineSm4(S);
XB(i)=S.y_hard(1);
LB(i)=S.x_hard(2)*2;
X(i)=XB(i)-max(S.y_mc);
Vry(1,i)=LB(i)/XB(i);
Vry(2,i)=X(i)/LB(i);


end
end
end
end
end
end


[~,idx] = sort(Vry(1,:)); % sort just the first column
Vry = Vry(:,idx); 
Vry(3,:)=0.62*(Vry(1,:).^-1.15);
Vry(4,:)=0.6784*(Vry(1,:).^-1.2148);

figure (14)
% plot(Vry(1,:),Vry(2,:))
% plot(Vry(1,:),Vry(2,:),'^',Vry(1,:),Vry(3,:),'linewidth','2','Color','c-',Vry(1,:),Vry(4,:))

plot(Vry(1,:),Vry(2,:),'^',Vry(1,:),Vry(3,:),Vry(1,:),Vry(4,:))
 xlabel('LB/XB');
 ylabel('X/LB');
 print('Analysis','-dpng', '-r300')
 legend('approach','KHOUNG','Silvester','Location',S.legendlocation);
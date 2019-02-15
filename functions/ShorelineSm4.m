function [S]=ShorelineSm(S0);
% MODEL : ShorelineS
%
% This model computes shoreline changes as a result of gradients in alongshore
% sediment transport for arbitrary shaped coastlines.
%
% INPUT:
%     S       data structure wiht fields:
%              .
%
% made by:      J.A. Roelvink (2016) - UNESCO-IHE
% modified by:  B.J.A. Huisman (2017) - Deltares

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFAULT INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S.Hso='';                                                                   % wave height [m]
S.phiw0=330*pi/180;                                                        % deep water wave angle [°N]
S.spread=90;                                                               % wave spreading [°] (wave_dir from range:  S.phiw0 +/- 0.5*S.spread)
S.WVCfile='';                                                              % wave time-series <-leave empty to use wave parameters ('S.Hso', 'S.phiw0' and 'S.spread')
S.Waveclimfile='';%'Waveclimate\wcon.txt';                                                         % wave climate file
S.Wavecorr=0;
S.d=10;                                                                    % active profile height [m]
S.ddeep=25;
S.dnearshore=8;
S.LDBcoastline='';%'linetest.xy';%Noord-Holland2007coastline.xy';                                                        % LDB with initial coastline shape ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
S.phif='';                                                         % Orientation of the foreshore [°N] <- only relevant for 'KAMP', 'MILH' or 'VR14'
S.trform='CERC';                                                           % switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
S.b=1e6;                                                                   % CERC : coeff in simple cerc formula
S.tper=6;                                                                  % KAMP & MILH : peak wave period [s]
S.d50=1.0e-3;                                                              % KAMP & MILH & VR14 : median grain diameter [m]
S.porosity=0.4;                                                            % KAMP & MILH & VR14 : S.porosity (typically 0.4) [-]
S.tanbeta=0.03;                                                            % KAMP & MILH & VR14 : mean bed slope [ratio 1/slope]
S.rhos=2650;                                                               % KAMP & MILH & VR14 : density of sand [kg/m3]
S.rhow=1025;                                                               % KAMP & MILH & VR14 : density of water [kg/m3]
S.g=9.81;                                                                  % KAMP & MILH & VR14 : gravitational acceleration [m2/s]
S.alpha=1.8;                                                               % KAMP & MILH & VR14 : calibration factor for point of breaking (S.alpha = 1.8 for Egmond data)
S.gamma=0.72;                                                              % KAMP & MILH & VR14 : breaking coefficient (Hs/h) with 5% breaking waves
S.Pswell=20;                                                               % VR14 : Percentage swell (between 0 - 100) [-]
S.crit=.9;                                                                 % stability criterion (not active)
S.reftime='';%'2017-01-01';                                                % Reference time (i.e. 'yyyy-mm-dd') <- leave empty to use t=0
%S.dt=.2;                                                                  % time step [year] (use 1/365/4 for time-series)
S.tc=1;
S.ds0=100;                                                                 % initial space step [m]
S.Courant=1;
S.ns=100;                                                                  % number of ... (not used)
%S.nt=5000;                                                                 % number of timesteps
S.twopoints=1;                                                             % switch for 'S.twopoints appraoch'
S.struct=0;                                                                % switch for using hard structures
S.LDBstructures='';                                                        % LDB with hard structures ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
S.growth=0;                                                                % calibration of nourishment growth rate
S.nourish=0;                                                               % switch (0/1) for nourishments
S.LDBnourish='';                                                           % LDB with nourishment locations ([Nx2] ASCII FILE WITHOUT HEADER) (i.e. polygon around relevant grid cells) <- leave empty to use interactive mode!
S.nourrate=100;
S.nourstart=0;
S.spit_width=50;                                                           % width of tip of spit
S.xlimits=[];                                                              % X limits of plot [m] [1x2] <- leave empty to automatically do this
S.ylimits=[];                                                              % Y limits of plot [m] [1x2] <- leave empty to automatically do this
S.XYwave =[];                                                              % X,Y location and scale of wave arrow [m] [1x2] (automatically determined on the basis of xy-limits and wave angle
S.XYoffset=[0,0];                                                          % shift in X,Y locaton for plotting <- leave empty to automatically shift grid <- use [0,0] for no shift
S.pauselength=[];                                                          % pause between subsequent timesteps (e.g. 0.0001) <- leave empty to not pause plot
S.outputdir='Output_Marina_del\';
S.rundir='';%'Delft3D\def_model\';
S.tide_interaction=false;
S.wave_interaction=false;
S.LDBplot = {};                                                            % cell array with at every line : string with LDB-filename, string with legend entry, string with plot format (e.g. 'b--') <- e.g. {'abc.ldb','line 1','k--'; 'def.ldb','line 2','r-.'; etc} <- leave empty to not use additional plots
S.legendlocation='northeast';
S.RWSfiletype=false;
S.fignryear=12;
S.plotinterval=1;
S.smoothfac=0;
S.endofsimulation='';%'2007-07-01';
ld=1000;                                                                    % Width of the land fill behind the shoreline [m]

%% boundary condition                                                       for (non cyclic) sections (ex.straight shoreline) ,If more than one (non cyclic)sections should adjusted manually
S.boundary_condition_start='FIXD2';                                               % boundary condition 'PRDC' for periodeic boundary condition or 'FIXD' for fixed boundary condition
S.boundary_condition_end='FIXD2';
S.BCfile='';                                                                % Boundary condition function in time (file) , columns (time , QS_start, Qs_end)
S.QS_start='';                                                              %[m3/year]
S.QS_end='';

%%Sources and Sinks
S.sources_sinks='';
S.SSfile='';

%% Extract shorelines
S.SLplot={};
S.extract_x_y=0;                                                            %get file with shorelines coordinates x,y
S.print_fig=0;
%% video
S.video=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET INPUT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fieldnms=fields(S0);
for ii=1:length(fieldnms)
    S.(fieldnms{ii}) = S0.(fieldnms{ii});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPARE COASTLINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x_mc,y_mc,x_mc0,y_mc0,S]=prepare_coastline(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPARE NOURISHMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if S.nourish
    if ~isempty(S.LDBnourish)
        xy_nour=load(S.LDBnourish);
        x_nour=xy_nour(:,1)'-S.XYoffset(1);
        y_nour=xy_nour(:,2)'-S.XYoffset(2);
    else
        [x_nour,y_nour]=select_multi_polygon('k');
    end
    [ x_n,y_n,n_nour,in1,in2 ] = get_one_polygon( x_nour,y_nour,1 );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPARE STRUCTURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x_hard,y_hard]=prepare_structures(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPARE WAVE CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[WVC,WC,timenum0]=prepare_wave_conditions(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAKE COMPUTATION (and plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all
figure(11);clf;
setFIGUREproperties(800,600,32);
n_mc=length(find(isnan(x_mc)))+1;
[x_mc,y_mc]=initialize_grid(x_mc,y_mc,S.ds0);
xlim(S.xlimits);ylim(S.ylimits);
%S.A=zeros(S.nt+1,10);
%iplot=0;
% for iw=1:size(WVC.timenum)
% indw(iw)=timenum0-WVC.timenum(iw)
% end
% indxw=find(min(indw)
tnow=timenum0;
tend = datenum(S.endofsimulation,'yyyy-mm-dd'); %HH:MM:SS

S.dt=0;
it=-1;   % for now = -1 to follow the exisiting code
vii=0;
%% For extracting specific shorelines
if ~isempty(S.SLplot)
    plot_time(:)=datenum(S.SLplot(:,1),'yyyy-mm-dd');
    tsl=1;
    tplot=plot_time(tsl);
    plot_time(end+1)=0;
else
    tplot=[];
    tsl=[];
    plot_time=[];
end
xmax=0;ymax=0;xmin=0;ymin=0;nmax=0;
xp_mc{:,:}={};
yp_mc{:,:}={};
vi = struct('cdata', cell(1,1), 'colormap', cell(1,1));



while  tnow<tend
    it=it+1;
    %% Introduce the wave conditions to the domain
    [phiw,S]=introduce_wave(S,WVC,WC,tnow);
    if S.wave_interaction
        [ Hg, phiwg ] = get_interpolated_wavefield_dir_Tp( S.xg,S.yg,S.Hg,S.dirg,S.Hso,phiw,S.tper,S.dirtab,S.Tptab);
    end
    
    %% Adaptive time step
%         [adt]=adaptive_time_step(x_mc,y_mc,S,phiw,Hg,phiwg);
    [adt]=adaptive_time_step(x_mc,y_mc,S,phiw);
    adt_record(it+1)=adt;
%     adt=625/(60*60*24*365);
    n_mc=length(find(isnan(x_mc)))+1;
    
    for i_mc=1:n_mc
        n_mc=length(find(isnan(x_mc)))+1;
        if i_mc>n_mc
            break
        end
        
        %% Alongshore coordinate s
        [s,x,y,x_mc,y_mc,ar]=make_sgrid_mc(x_mc,y_mc,S.ds0,i_mc,S.smoothfac);
        if S.wave_interaction
            [H phiw_cd]=get_refracted_waves(x,y,500,S.xg,S.yg,Hg,phiwg);
            phiwi=phiw_cd;
            S.Hso=H;
        else
            phiwi=phiw;
        end
        
        %pause
        itt=it+1;
        S.A(itt,i_mc)=double(ar);
        n=length(x)-1;
        %% Cyclic or not ?
        cyclic = hypot(x(end)-x(1),y(end)-y(1))<S.ds0;
        %% Angles
        [philoc,thetacrit,hsbr,sphibr,hbr]=angles(S,s,x,y,n,phiwi);
        %% Shadowing effect on Transport
        [xS,yS,shadowS,shadowS_h,shadow,philoc_cor,sphibr,S,hsbr]=transport_shadow_treat(x,y,x_mc,y_mc,x_hard,y_hard,phiw,philoc,sphibr,S,hsbr);
        %% Long shore Transport
        [QS,QSmax]=transport(S,it,s,philoc_cor,hsbr,sphibr);
        %% wave diffraction
        if ~cyclic
     %       [QS,xp,yp]=wave_diffraction11(QS,x,y,S,hsbr,sphibr,phiw,x_mc,y_mc,x_hard,y_hard,n,hbr);
        end
        %% Upwind correction for high-angle
        [QS]=upwind_correction(philoc,n,cyclic,S,thetacrit,shadowS,shadowS_h,QSmax,QS);
        
        %% Set QS=0 for large angles (>90°)
        QS(abs(philoc_cor)>pi/2)=0.;
        
        %               if length(shadowS_h)>0 % to be tested with structure case
        %                  QS(shadowS_h)=0;
        %             end
        %          QS(shadowS)=0;
        
        %% Coastline cells intersected by hard structures
        if length(x_hard)>1
            [structS]=find_struct(x,y,x_hard,y_hard);
            QS(structS)=0;
        else
            structS=[];
        end
        %% Nourisment?
        nour=zeros(size(x));
        if S.nourish
            if it>S.nourstart
                for i_n=1:n_nour
                    [ x_n,y_n,n_nour,in1,in2 ] = get_one_polygon( x_nour,y_nour,i_n );
                    nour=nour+inpolygon(x,y,x_n,y_n);
                end
            end
        end
        
        
        %% boundary condition
        [QS]=Boundary_condition(QS,S,cyclic,tnow,it);
        %% Sources and Sinks
        [QS]=Sources_Sinks(QS,S,tnow,x_mc,y_mc);
        
        %% Coastline change
        [x,y,S]=coastline_change(S,s,QS,x,y,cyclic,n,tplot,tnow,tend,adt,nour,it,WVC);
        
        [ x,y, spit,width] = find_spit_width_v2( s,x,y,phiw,shadow,S.spit_width );
        %disp(['max(dx)',num2str(max(dx)),' max(x)',num2str(max(x))])
        %% insert x and y back into x_mc,y_mc
        [x_mc,y_mc]=insert_section(x,y,x_mc,y_mc,i_mc);
        %% Merge coastlines where necessary
        [ xnew,ynew ]=merge_coastlines( x,y );
        [x_mc,y_mc]=insert_section(xnew,ynew,x_mc,y_mc,i_mc);
        %% clean up redundant NaNs
        [x_mc,y_mc]=cleanup_nans(x_mc,y_mc);
        if 0
            plot(x_mc0,y_mc0,x_mc,y_mc,'.r',xS(shadowS),yS(shadowS),'.g',xS(structS),yS(structS),'*k',x_hard,y_hard,'k', 'linewidth',2);
            xlabel('Easting [m]');
            ylabel('Northing [m]');
            %plot(x_mc0,y_mc0,x_mc,y_mc,'.-r',x_hard,y_hard,'k', 'linewidth',2)
            axis equal;
            hold on;
            plot(xS(shadowS_h),yS(shadowS_h),'bo');
            if 1
                arx=[-2000,-2000+500*cos(3/2*pi-phiw)];
                ary=[4000,4000+500*sin(3/2*pi-phiw)];
                plot(arx,ary,'linewidth',2);
            end
            title([num2str(it)]);
            hold off;
            xlim(S.xlimits);
            ylim(S.ylimits);
            drawnow;
            fname=['cl',num2str(1000+it),'png'];
            %if mod(it,10)==0
            %    print('-dpng',fname)
            %end
            if ~isempty(S.pauselength)
                pause(S.pauselength);
            end
        end
    end
    [ x_mc,y_mc ] = merge_coastlines_mc(x_mc,y_mc);
    %% clean up redundant NaNs
    [x_mc,y_mc]=cleanup_nans(x_mc,y_mc);
    %% plot settings
    
%     [vi,xp_mc,yp_mc,S,xmax,ymax,xmin,ymin,nmax,tsl,tplot,vii]=PlotX(S,x_hard,y_hard,x_mc0,y_mc0,it,x_mc,y_mc,ld,tnow,phiw,tplot,tsl,plot_time,xmax,ymax,xmin,ymin,nmax,xp_mc,yp_mc,vi,vii,xp,yp,QS,xS);
            [vi,xp_mc,yp_mc,S,xmax,ymax,xmin,ymin,nmax,tsl,tplot,vii]=PlotX(S,x_hard,y_hard,x_mc0,y_mc0,it,x_mc,y_mc,ld,tnow,phiw,tplot,tsl,plot_time,xmax,ymax,xmin,ymin,nmax,xp_mc,yp_mc,vi,vii);
    
    S.x_mc=x_mc;
    S.y_mc=y_mc;
    S.x_mc0=x_mc0;
    S.y_mc0=y_mc0;
    % S.timenum0=timenuma(end);
    
    % Apply shoreline change due to tide
    if S.tide_interaction
        [S]=update_shoreline(S);
    end
    tnow = tnow+ S.dt*365; %dt [year]   %calculate the current time after each time step
    S.times(it+1)=tnow;
end

%% Extract shorelines cooridnates /figures at specific dates
extract_shoreline(S,x_hard,y_hard,x_mc0,y_mc0,xp_mc,yp_mc);
%% videos
make_video(S,vi);

% print((strcat('xhrd_',num2str(S.x_hard),'_yhrd_',num2str(S.y_hard),'2')),'-dpng', '-r300')

% figure(13)
% xit=[1:1:length(adt_record)];
% plot(xit,adt_record,'LineWidth',2)
% xlabel('Steps');
%  ylabel('Adaptive time steps (years)');
%  print('adt','-dpng', '-r300')
%  save('adt2','adt_record','-ascii');
%  figure(14)
% 
% plot(xit,adt_record,'*')
% xlabel('Steps');
%  ylabel('Adaptive time steps (second)');
%  print('adt2','-dpng', '-r300')


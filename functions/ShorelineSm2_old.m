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
S.Hso=1;                                                                   % wave height [m]
S.phiw0=330*pi/180;                                                        % deep water wave angle [°N]
S.spread=90;                                                               % wave spreading [°] (wave_dir from range:  S.phiw0 +/- 0.5*S.spread)
S.WVCfile='';                                                              % wave time-series <-leave empty to use wave parameters ('S.Hso', 'S.phiw0' and 'S.spread')
S.Waveclimfile='Waveclimate\wcon.txt';                                                         % wave climate file
S.Wavecorr=0;
S.d=10;                                                                    % active profile height [m]
S.ddeep=25;
S.dnearshore=8;
S.LDBcoastline='Noord-Holland2007coastline.xy';                                                        % LDB with initial coastline shape ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
S.phif=330*pi/180;                                                         % Orientation of the foreshore [°N] <- only relevant for 'KAMP', 'MILH' or 'VR14'
S.trform='CERC';                                                           % switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
S.b=1e6;                                                                   % CERC : coeff in simple cerc formula
S.tper=6;                                                                  % KAMP & MILH : peak wave period [s]
S.d50=2.0e-4;                                                              % KAMP & MILH & VR14 : median grain diameter [m]
S.porosity=0.4;                                                            % KAMP & MILH & VR14 : S.porosity (typically 0.4) [-]
S.tanbeta=0.03;                                                            % KAMP & MILH & VR14 : mean bed slope [ratio 1/slope]
S.rhos=2650;                                                               % KAMP & MILH & VR14 : density of sand [kg/m3]
S.rhow=1025;                                                               % KAMP & MILH & VR14 : density of water [kg/m3]
S.g=9.81;                                                                  % KAMP & MILH & VR14 : gravitational acceleration [m2/s]
S.alpha=1.8;                                                               % KAMP & MILH & VR14 : calibration factor for point of breaking (S.alpha = 1.8 for Egmond data)
S.gamma=0.72;                                                              % KAMP & MILH & VR14 : breaking coefficient (Hs/h) with 5% breaking waves
S.Pswell=20;                                                               % VR14 : Percentage swell (between 0 - 100) [-]
S.crit=.9;                                                                 % stability criterion (not active)
S.reftime='2017-01-01';                                                               % Reference time (i.e. 'yyyy-mm-dd') <- leave empty to use t=0
%S.dt=.2;                                                                  % time step [year] (use 1/365/4 for time-series)
S.ds0=100;                                                                 % initial space step [m]
S.Courant=1;
S.ns=100;                                                                  % number of ... (not used)
%S.nt=5000;                                                                 % number of timesteps
S.twopoints=1;                                                             % switch for 'S.twopoints appraoch'
S.struct=1;                                                                % switch for using hard structures
S.LDBstructures='alt2.xy';                                                        % LDB with hard structures ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
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
S.endofsimulation='2027-01-01';
ld=3000;                                                                    % Width of the land fill behind the shoreline [m]

%% boundary condition                                                       for (non cyclic) sections (ex.straight shoreline) ,If more than one (non cyclic)sections should adjusted manually 
S.boundary_condition_start='';                                               % boundary condition 'PRDC' for periodeic boundary condition or 'FIXD' for fixed boundary condition
S.boundary_condition_end='';
S.BCfile='';                                                                % Boundary condition function in time (file) , columns (time , QS_start, Qs_end)
S.QS_start='';                                                              %[m3/year]
S.QS_end='';

%%Sources and Sinks
S.sources_sinks='';
S.SSfile='';

%% Extract shorelines
S.SLplot={};



% cell array with at every line : string with LDB-filename, string with legend entry, string with plot format (e.g. 'b--') <- e.g. {'abc.ldb','line 1','k--'; 'def.ldb','line 2','r-.'; etc} <- leave empty to not use additional plots


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
if ~isempty(S.LDBcoastline)&~isfield(S,'x_mc')
    xy_mc=load(S.LDBcoastline);
    if isempty(S.XYoffset)
        S.XYoffset = [floor(min(xy_mc(:,1))/1000)*1000 , floor(min(xy_mc(:,2))/1000)*1000];
    end
    x_mc=xy_mc(:,1)' - S.XYoffset(1);   %x_mc=xy_mc(end:-1:1,1)';   % SHIFT COASLTINE
    y_mc=xy_mc(:,2)' - S.XYoffset(2);   %y_mc=xy_mc(end:-1:1,2)';   % SHIFT COASLTINE
    x_mc0=x_mc;
    y_mc0=y_mc;
    
elseif ~isfield(S,'x_mc')
    figure(99);clf;
    plot([S.xlimits(1) S.xlimits(2) S.xlimits(2) S.xlimits(1) S.xlimits(1)], ...
        [S.ylimits(1) S.ylimits(1) S.ylimits(2) S.ylimits(2) S.ylimits(1)],'k:');
    axis equal;
    xl=xlim;yl=ylim;
    xlabel('Easting [m]');
    ylabel('Northing [m]');
    htxt=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add coastline (LMB); Next segment (RMB); Exit (q)');set(htxt,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.1 0.6]);
    [x_mc,y_mc]=select_multi_polygon('r');
    set(htxt,'Visible','off');
    x_mc0=x_mc;
    y_mc0=y_mc;
    
else
    x_mc=S.x_mc;
    y_mc=S.y_mc;
    x_mc0=S.x_mc0;
    y_mc0=S.y_mc0;
    
end
if isempty(S.xlimits)
    S.xlimits = [min(x_mc),max(x_mc)];
else
    S.xlimits = S.xlimits - S.XYoffset(1);
end
if isempty(S.ylimits)
    S.ylimits = [min(y_mc),max(y_mc)];
else
    S.ylimits = S.ylimits - S.XYoffset(2);
end

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
if S.struct
    if isfield(S,'x_hard')
        x_hard=S.x_hard;
        y_hard=S.y_hard;
    elseif ~isempty(S.LDBstructures)
        xy_hard=load(S.LDBstructures);
        x_hard=xy_hard(:,1)'-S.XYoffset(1);
        y_hard=xy_hard(:,2)'-S.XYoffset(2);
        S.x_hard=x_hard;
        S.y_hard=y_hard;
    else
        figure(99);
        axis equal;
        xl=xlim;yl=ylim;
        htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add structure (LMB); Next structure (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
        [x_hard,y_hard]=select_multi_polygon('k');
        set(htxt2,'Visible','off');
        S.x_hard=x_hard;
        S.y_hard=y_hard;
    end
else
    x_hard=[];
    y_hard=[];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PREPARE WAVE CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(S,'timenum0')
    if ~isempty(S.reftime)
        timenum0 = datenum(S.reftime,'yyyy-mm-dd'); %HH:MM:SS
    else
        timenum0 = 0;
    end
else
    timenum0=S.timenum0;
end
%timenum = timenum0+[0:S.nt]*S.dt*365;
if ~isempty(S.WVCfile)
    if S.RWSfiletype
        WVCraw=load(S.WVCfile);warning off;
        WVC=struct;
        WVC.timenum=datenum([num2str(WVCraw(:,1))],'yyyymmddHHMM'); %,'%08.0f'),num2str(WVCraw(:,2),'%06.0f')
        WVC.Hs=interpNANs(WVCraw(:,2));
        WVC.Tp=interpNANs(WVCraw(:,3));
        WVC.Dir=interpNANsDIR(WVCraw(:,4));
    else
        WVCraw=load(S.WVCfile);warning off;
        WVC=struct;
        WVC.timenum=datenum([num2str(WVCraw(:,1),'%08.0f'),num2str(WVCraw(:,2),'%06.0f')],'yyyymmddHHMMSS');
        WVC.Hs=interpNANs(WVCraw(:,2));
        WVC.Tp=interpNANs(WVCraw(:,3));
        WVC.Dir=interpNANsDIR(WVCraw(:,5));
    end
timenuma=WVC.timenum;
%     WVC.timenumi=timenum;
%     WVC.Hsi=interp1(WVC.timenum,WVC.Hs,timenum);
%     WVC.Tpi=interp1(WVC.timenum,WVC.Tp,timenum);
%     WVC.Diri=mod(atan2(interp1(WVC.timenum,sin(WVC.Dir*pi/180),timenum),interp1(WVC.timenum,cos(WVC.Dir*pi/180),timenum))*180/pi,360);
%     WVC.Hsi=interpNANs(WVC.Hsi);
%     WVC.Tpi=interpNANs(WVC.Tpi);
%     WVC.Diri=interpNANsDIR(WVC.Diri);
%     figure;plot(WVC.timenum,WVC.Dir,'b.',WVC.timenumi,WVC.Diri,'k.');
%     pause

elseif ~isempty(S.Waveclimfile)
    WCraw=load(S.Waveclimfile);
    WC=struct;
    WC.Hs=WCraw(:,1);
    WC.Tp=WCraw(:,2);
    WC.dir=WCraw(:,3)+S.Wavecorr;
    timenuma=timenum0;
else
     timenuma=timenum0;  
end


%%%%%
%% 


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
iplot=0;

tnow=timenum0;       
timenume = datenum(S.endofsimulation,'yyyy-mm-dd'); %HH:MM:SS

S.dt=0;
it=-1;   % for now = -1 to follow the exisiting code   

%% For extracting specific shorelines
if ~isempty(S.SLplot)
    plot_time(:)=datenum(S.SLplot(:,1),'yyyy-mm-dd');
    tsl=1;
    tplot=plot_time(tsl);
    plot_time(end+1)=0;
end

while  tnow<timenume


    % tnow = timenum0 + (it)*S.dt*365;
    it=it+1;
    %tnow = tnow+ S.dt*365; %dt [year]   %calculate the current time after each time step




%the following lines removed from (198 :206) to inside the loop
if ~isempty(S.WVCfile)

       %WVC.timenumi=timenuma;  %follow the consequences 
    WVC.Hsi=interp1(WVC.timenum,WVC.Hs,tnow); %If the interpolation needs more computational .. We could interpolate in between 
    WVC.Tpi=interp1(WVC.timenum,WVC.Tp,tnow);
    WVC.Diri=mod(atan2(interp1(WVC.timenum,sin(WVC.Dir*pi/180),tnow),interp1(WVC.timenum,cos(WVC.Dir*pi/180),tnow))*180/pi,360);
    S.Hso=interpNANs(WVC.Hsi);
    S.tper=interpNANs(WVC.Tpi);
    phiw=interpNANsDIR(WVC.Diri);
    %figure;plot(WVC.timenum,WVC.Dir,'b.',WVC.timenumi,WVC.Diri,'k.');
   % pause
end
  
n_mc=length(find(isnan(x_mc)))+1;
    for i_mc=1:n_mc
        if ~isempty(S.WVCfile)
            phiw=WVC.Diri*pi/180;
            S.Hso=WVC.Hsi;
            S.tper=WVC.Tpi; 
        elseif ~isempty(S.Waveclimfile)
            iwc=round((rand*length(WC.Hs)+0.5)); 
            phiw=WC.dir(iwc)*pi/180;
            S.Hso=WC.Hs(iwc);
            S.tper=WC.Tp(iwc);
        else
            phiw=S.phiw0+(rand-.5)*pi/180*S.spread;
            S.Hs0=1;
        end
        if S.wave_interaction
            [ Hg, phiwg ] = get_interpolated_wavefield( S.xg,S.yg,S.Hg_all,S.phiwg_all,S.Hso,phiw*180/pi);
        end
        
        n_mc=length(find(isnan(x_mc)))+1;
        if i_mc>n_mc
            break
        end
        %% Alongshore coordinate s
        [s,x,y,x_mc,y_mc,ar]=make_sgrid_mc(x_mc,y_mc,S.ds0,i_mc,S.smoothfac);
        if S.wave_interaction
            [H phiw_cd]=get_refracted_waves(x,y,500,S.xg,S.yg,Hg,phiwg);
        end
        %pause
        itt=it+1;
        S.A(itt,i_mc)=double(ar);
        n=length(x)-1;
        %% Cyclic or not ?
        cyclic = hypot(x(end)-x(1),y(end)-y(1))<S.ds0;
        %% Angles
        sphi=[];
        phic=[];
        for i=1:n
            sphi(i)=.5*(s(i)+s(i+1));
            phic(i)=2*pi-atan2(y(i+1)-y(i),x(i+1)-x(i));
        end

        if S.wave_interaction
            philoc=atan2(sin(phic-phiw_cd),cos(phic-phiw_cd));
        else
            philoc=atan2(sin(phic-phiw),cos(phic-phiw));
        end
        if abs(mean(philoc))>pi/2
            continue
        end
        thetacrit0 = repmat(pi/4,[1 n]);
        thetacrit=thetacrit0;
        
        %% Transform to nearshore
        if strcmpi(S.trform,'CERC')
            thetacrit=thetacrit0;
        elseif strcmpi(S.trform,'CERC2')
            %           phi9=[0:.01:70]*pi/180;id=find(cos(phi9).^(6/5).*sin(phi9)==max(cos(phi9).^(6/5).*sin(phi9)));
            %           thetacrit=thetacrit0*(phi9(id)*180/pi/45);
            [kh0,c0]=GUO2002(S.tper,S.ddeep);
            [kh1,c1]=GUO2002(S.tper,S.dnearshore);
            thetacrit=asin((c1./c0).*sin(thetacrit0));
        else
            % refraction on foreshore (S.ddeep -> S.dnearshore)
            if isempty(S.phif)
                S.phif=phic;
            end
            philoc0=atan2(sin(S.phif-phiw),cos(S.phif-phiw));
            [kh0,c0]=GUO2002(S.tper,S.ddeep);
            [kh1,c1]=GUO2002(S.tper,S.dnearshore);
            philoc1=asin((c1./c0).*sin(philoc0));
            thetacrit1=asin((c1./c0).*sin(thetacrit0));
            
            % refraction in the nearshore
            philoc1B = philoc1+(phic-S.phif);
            hbr=real(((S.Hso.^2 .* c1 .* cos(philoc))./(S.alpha .* S.gamma.^2 .* S.g.^0.5)).^0.4);
            hsbr=hbr*S.gamma;
            [kbr,cbr]=GUO2002(S.tper,hbr);
            philoc2=asin((cbr./c1).*sin(philoc1B));
            thetacrit2=asin((cbr./c1).*sin(thetacrit1));
            
            thetacrit=thetacrit2;
            philoc=philoc2;
        end
        
        
        
        
        %% Upwind correction for high-angle
        philoc_cor=philoc;
        %[ xS,yS,shadowS ]        = find_shadows_mc( x,y,[x_mc,nan,x_hard],[y_mc,nan,y_hard],phiw );
        [ xS,yS,shadowS ]        = find_shadows_mc( x,y,[x_mc],[y_mc],phiw,0 );
        if length(x_hard)>0&length(x)>0
            [ xS,yS,shadowS_h ]      = find_shadows_mc( x,y,[x_hard],[y_hard],phiw,1 );
            shadow=zeros(size(x));
            shadow(1:end-1)=shadowS_h;
            shadow(2:end)=min(shadow(2:end),shadowS_h);
        else
            shadowS_h=[];
            shadow=zeros(size(x));
        end
                philoc_cor(shadowS)=0;
               philoc_cor(shadowS_h)=0;
        
        for i=1:n
            if cyclic
                im1=mod2(i-1,n);
                im2=mod2(i-2,n);
                ip1=mod2(i+1,n);
                ip2=mod2(i+2,n);
            else
                im1=max(i-1,1);
                im2=max(i-2,1);
                ip1=min(i+1,n);
                ip2=min(i+2,n);
            end
            %if abs(philoc(i))>pi/4&&philoc(im1)<pi/4&&philoc(im1)>0
            if philoc(i)>thetacrit(i)&&philoc(im1)<thetacrit(i)&& ...
                    philoc(im1)>0&&~shadowS(im1)&& ...
                    (isempty(shadowS_h)||(~isempty(shadowS_h)&&~shadowS_h(ip1)&&~shadowS_h(ip2)))
                philoc_cor(i)=thetacrit(i);
                if S.twopoints
                    %philoc_cor(ip1)=pi/15;
                    philoc_cor(ip1)=thetacrit(i)/2;
                    philoc_cor(ip2)=0;
                else
                    philoc_cor(ip1)=0;
                end
                %elseif abs(philoc(i))>pi/4&&philoc(ip1)>-pi/4&&philoc(ip1)<0
            elseif philoc(i)<-thetacrit(i)&&philoc(ip1)>-thetacrit(i)&& ...
                    philoc(ip1)<0&&~shadowS(ip1)&& ...
                    (isempty(shadowS_h)||(~isempty(shadowS_h)&&~shadowS_h(im1)&&~shadowS_h(im2)))
                philoc_cor(i)=-thetacrit(i);
                if  S.twopoints
                    %philoc_cor(im1)=-pi/15;
                    philoc_cor(im1)=-thetacrit(i)/2;
                    philoc_cor(im2)=0;
                else
                    philoc_cor(im1)=0;
                end
            end
        end
        
       %  philoc_cor(1)=OSd;
      %  philoc_cor(end)=OSd2;
     
      if it==0   %four boundary conditions
          QS=zeros(size(s));
          S.QS_start=QS(1);
          S.QS_end=QS(end);
      end
        %% Transport : CERC
        if strcmpi(S.trform,'CERC')
             k=0.2;                                                      % using CERC (1984) value of k(SPM,Hs)=0.39
            % b_theory = k * (S.rhow * S.g^0.5 / (16 * sqrt(S.gamma)* (S.rhos-S.rhow) * (1-S.porosity)));  % b_theory = 0.0946 * 365*24*60*60 = 2.9833E+6
            QS=S.b*S.Hso^2.5*sin(2*(philoc_cor));QScerc=QS;
            %Smax=S.b*S.Hso^2.5;
        %adt=0.25*S.ds0^2/(S.b*S.Hso^2.5*cos(2*(phiw)))*S.d%/(365*24*60*60)^2 %should it be S.ds0 ?
         adt=0.25*S.ds0^2/(S.b*S.Hso^2.5*1)*S.d %(365*24*60*60)^2 %assumption cos2deltb=1
        %adt=(0.25*S.ds0^2)/(k * (S.rhow * S.g^0.5 / (16 * sqrt(S.gamma)* (S.rhos-S.rhow) * (1-S.porosity))))*S.d/(365*24*60*60)
        end
        if strcmpi(S.trform,'CERC2')
            k=0.12;                                                      % using CERC (1984) value of k(SPM,Hs)=0.39
            b1 = k * (S.rhow * S.g^0.5 / (16 * sqrt(S.gamma)* (S.rhos-S.rhow) * (1-S.porosity)));  % b_theory = 0.0946 * 365*24*60*60 = 2.9833E+6
            b2 = b1 * ((S.gamma.*S.g).^0.5 /(2*pi)).^0.2;  % b_theory = 0.0946 * 365*24*60*60 = 2.9833E+6
            QS=365*24*60*60*b2*S.Hso.^(12/5)*S.tper.^(1/5).*cos(philoc_cor).^(6/5).*sin(philoc_cor);QScerc2=QS;
            %Smax=S.b*S.Hso^2.5;
            %adt=0.5*S.ds0^2/(-b2/S.d*S.Hso.^(12/5)*S.tper.^(1/5).*cos(philoc_cor).^(1/5)*((cos(philoc_cor).^2)-(6/5)*(sin(philoc_cor).^2)))/(365*24*60*60);
             adt=(0.5*(S.ds0.^2))./((b2/S.d*S.Hso.^(12/5)).*S.tper.^(1/5))/(365*24*60*60);%*cos(philoc_cor).^(1/5).*((cos(philoc_cor).^2)-(6/5).*(sin(philoc_cor).^2)))/(365*24*60*60);
            % adt=(0.5*(S.ds0.^2))./((b2/S.d*S.Hso.^(12/5)).*S.tper.^(1/5)*cos(philoc_cor).^(1/5).*((cos(philoc_cor).^2)+(6/5).*(sin(philoc_cor).^2)))/(365*24*60*60);
            % adt=min(adt)
             %adt=min(adt(adt > 0))
            % adt=real(adt)
            
        end
        if ~strcmpi(S.trform,'CERC') && ~strcmpi(S.trform,'CERC2')
            %[kh0,c0]=GUO2002(S.tper,S.d);
            %hbr=real(((S.Hso.^2 .* c0 .* cos(philoc))./(S.alpha .* S.gamma.^2 .* S.g.^0.5)).^0.4);
            %[kbr,cbr]=GUO2002(S.tper,hbr);
            %sphibr=asin((cbr./c0).*sin(philoc));
            %hsbr=hbr*S.gamma;
            
            % refraction on foreshore (S.ddeep -> S.dnearshore)
            if isempty(S.phif)
                S.phif=phic;
            end
            philoc0=atan2(sin(S.phif-phiw),cos(S.phif-phiw));
            [kh0,c0]=GUO2002(S.tper,S.ddeep);
            [kh1,c1]=GUO2002(S.tper,S.dnearshore);
            philoc1=asin((c1./c0).*sin(philoc0));
            
            % refraction in the nearshore
            philoc1B = philoc1+(philoc_cor-philoc0);%(phic-S.phif);
            hbr=real(((S.Hso.^2 .* c1 .* cos(philoc))./(S.alpha .* S.gamma.^2 .* S.g.^0.5)).^0.4);
            hsbr=hbr*S.gamma;
            [kbr,cbr]=GUO2002(S.tper,hbr);
            sphibr=asin((cbr./c1).*sin(philoc1B));
        end
        %% Transport : Kamphuis
        if strcmpi(S.trform,'KAMP')
            QSkampmass=2.33 * S.rhos/(S.rhos-S.rhow) .* S.tper.^1.5 .* S.tanbeta.^0.75 .* S.d50.^-0.25 .* hsbr.^2 .* real((sin(2*sphibr)).^0.6);
            QS = 365*24*60*60*(QSkampmass /(S.rhos-S.rhow)) /(1.0-S.porosity);QSkamp=QS;
           % adt=0.02
           adt=(0.5*(S.ds0.^2))./(QS).*S.d.*sphibr;
           adt=min(adt(adt > 0))
            %% Transport : Mil-Homens
        elseif strcmpi(S.trform,'MILH')
            QSmilhmass=0.15*S.rhos/(S.rhos-S.rhow) .* S.tper.^0.89 .* S.tanbeta.^0.86 .* S.d50.^-0.69 .* hsbr.^2.75 .* real((sin(2*sphibr)).^0.5);
            QS = 365*24*60*60*(QSmilhmass /(S.rhos-S.rhow)) /(1.0-S.porosity);QSmilh=QS;
            %% Transport : Van Rijn (2014)
        elseif strcmpi(S.trform,'VR14')
            kswell=0.015*S.Pswell+(1-0.01*S.Pswell);
            vwave=0.3*real(sin(2*sphibr)).*(S.g.*hsbr).^0.5;
            vtide=0;
            vtotal=vwave+vtide;
            QSvr14mass=0.0006 .* kswell .* S.rhos .* S.tanbeta.^0.4 .*S.d50.^-0.6 .* hsbr.^2.6 .* abs(vtotal);
            QS = 365*24*60*60*(QSvr14mass /(S.rhos-S.rhow)) /(1.0-S.porosity);QSvr14=QS;
            %%Transport: CERC (vitousek,Barnard2015)
        elseif strcmpi(S.trform,'CERC3')
            QS=S.b*hsbr.^2.5.*sin(2*(philoc_cor));QScerc3=QS;
            adt=0.25*S.ds0.^2./(S.b*hsbr.^2.5)*S.d;
            adt=min(adt)
        end
        
        
        %% Set QS=0 for large angles (>90°)
        QS(abs(philoc_cor)>pi/2)=0.;
        if length(shadowS_h)>0
            QS(shadowS_h)=0;
         end
        %% Coastline cells intersected by hard structures
        if length(x_hard)>0
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
        dSds=zeros(size(s));
        ndev=zeros(size(s));
        dndt=zeros(size(s));
        dx=zeros(size(s));
        dy=zeros(size(s));
        for i=1:n
            if cyclic
                im1=mod2(i-1,n);
                ip1=mod2(i+1,n);
            else
                im1=max(i-1,1);
                ip1=min(i+1,n+1);
            end
            
            
            dSds(i)=(QS(i)-QS(im1))*2/max(hypot(x(ip1)-x(im1),y(ip1)-y(im1)),S.ds0);
            %dSds(i)=(QS(i)-QS(im1))*2/hypot(x(ip1)-x(im1),y(ip1)-y(im1));
            %dSds(i)=(S(i)-S(im1))*2/...
            %   (hypot(x(ip1)-x(i),y(ip1)-y(i))+hypot(x(i)-x(im1),y(i)-y(im1)));
            dQ(i)=(QS(i)-QS(im1));
           
            
%            
           
            dndt(i)=-dSds(i)/S.d;
            ndev(i)=hypot(x(i)-.5*(x(im1)+x(ip1)),y(i)-.5*(y(im1)+y(ip1)));
        end
        
        
        % Bondary Conditions
           if strcmpi(S.boundary_condition_start,'PRDC') && ~cyclic
               dSds(1)=(QS(1)-QS(end))*2/max((hypot(x(2)-x(1),y(2)-y(1))+hypot(x(end)-x(end-1),y(end)-y(end-1))),S.ds0);    
           elseif strcmpi(S.boundary_condition_start,'PRDC2') && ~cyclic
               dndt(1)=dndt(end-1);
               elseif strcmpi(S.boundary_condition_start,'FIXD2') && ~cyclic
    dndst(1)=dndt(2);
           end
           if strcmpi(S.boundary_condition_end,'PRDC') && ~cyclic
               dSds(end)=(QS(end)-QS(end-1))*2/max((hypot(x(2)-x(1),y(2)-y(1))+hypot(x(end)-x(end-1),y(end)-y(end-1))),S.ds0);
           elseif strcmpi(S.boundary_condition_end,'PRDC2') && ~cyclic
               dndt(end)=dndt(2);
               elseif strcmpi(S.boundary_condition_end,'FIXD2') && ~cyclic
    dndt(end)=dndt(end-1);
           end
        for i=1:n
            if cyclic
                im1=mod2(i-1,n);
                ip1=mod2(i+1,n);
            else
                im1=max(i-1,1);
                ip1=min(i+1,n+1);
            end

            if ~isempty(S.SLplot) && tplot~=0 %%for shorelines extraction
            
               if ~isempty(S.WVCfile)
                S.dt=min(S.tc*adt,(min((WVC.timenum(it+2)-tnow)/365,min((timenume-tnow)/365,(tplot-tnow)/365))));
            else
                S.dt=min(S.tc*adt,min((timenume-tnow)/365,(tplot-tnow)/365));
            end 
                
            else
                if ~isempty(S.WVCfile)
                S.dt=min(S.tc*adt,(min((WVC.timenum(it+2)-tnow)/365,(timenume-tnow)/365)));
            else
                S.dt=min(S.tc*adt,(timenume-tnow)/365);
            end
            end




            dn=(dndt(i)+S.growth+nour(i)*S.nourrate)*S.dt;
            
           
           
            dx(i)=-dn*(y(ip1)-y(im1))/hypot(x(ip1)-x(im1),y(ip1)-y(im1));
            dy(i)= dn*(x(ip1)-x(im1))/hypot(x(ip1)-x(im1),y(ip1)-y(im1));
            %courant=S.dt/S.ds0/S.d*dQ/dn
             adt2(i)=S.ds0*S.d*dn/dQ(i)*S.Courant;
        end
         adtc=abs(min(adt2));
         
          if abs(min(adt2))*1.00000001 < S.dt
                'HEEY'
               %pause
            end
        for i=1:n
            x(i)=x(i)+dx(i);
            y(i)=y(i)+dy(i);    
        end
        if cyclic
            x(n+1)=x(1);
            y(n+1)=y(1);
        end
        
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
    if mod(it,S.plotinterval)==0
        iplot=iplot+1;
        if 1
            %subplot(121)
            if it == 0 %for plottng purpose
            xmax=0;ymax=0;xmin=0;ymin=0;nmax=0;
            end
            [xmax,ymax,xmin,ymin,nmax]=fill_sections(x_mc,y_mc,ld,it,xmax,ymax,xmin,ymin,nmax);
            %axis([-10000 10000 0 10000])
            %set(gca,'xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[])
            hold on;
            plot(x_hard,y_hard,'k','linewidth',2);
            plot(x_mc0,y_mc0,'b','linewidth',2);
            if isempty(S.XYwave)
                S.XYwave = [S.xlimits(1)+(sin(S.phiw0)+1)/2*diff(S.xlimits),S.ylimits(1)+(cos(S.phiw0)+1)/2*diff(S.ylimits),diff(S.xlimits)/12];
            end
            if ~isempty(S.reftime) 
                text(S.XYwave(1)-S.XYoffset(1),S.XYwave(2)-S.XYoffset(2)+diff(S.ylimits)/12,['H_{m0}=',num2str(S.Hso,'%2.1f'),' (',datestr(tnow,'yyyy-mmm'),')']);
            else
                text(S.XYwave(1)-S.XYoffset(1),S.XYwave(2)-S.XYoffset(2)+diff(S.ylimits)/12,['H_{m0}=',num2str(S.Hso,'%2.1f'),' (',num2str(tnow/365,'%3.2f'),' yr)']);
            end
            if 1
                arx=[S.XYwave(1)-S.XYoffset(1),S.XYwave(1)-S.XYoffset(1)+S.XYwave(3)*cos(3/2*pi-phiw)];
                ary=[S.XYwave(2)-S.XYoffset(2),S.XYwave(2)-S.XYoffset(2)+S.XYwave(3)*sin(3/2*pi-phiw)];
                %plot(arx,ary,'linewidth',2);
                qvrscale=1.5;
                hqvr=quiver(arx(1),ary(1),S.XYwave(3)*cos(3/2*pi-phiw),S.XYwave(3)*sin(3/2*pi-phiw),qvrscale);
                set(hqvr,'linewidth',2,'Color','k','AutoScale','on','AutoSCaleFactor',1.5);
                %set(hqvr,'MaxHeadSize',50,'MarkerSize',100);
            end
            
            %% PLOT REFERENCE LINES AND LEGEND
            if ~isempty(S.LDBplot)
                for mm=1:size(S.LDBplot,1)
                    LDBplotval = landboundary('read',S.LDBplot{mm,1});
                    hp6(mm)=plot(LDBplotval(:,1)-S.XYoffset(1),LDBplotval(:,2)-S.XYoffset(2),S.LDBplot{mm,3},'linewidth',1.5);
                end
                %hleg = legend(hp6,S.LDBplot(:,2)','Location',S.legendlocation);
               % set(hleg,'Box','off','Color','None');
            end
            
            %% Plot shorelines at specific dates
            
             if ~isempty(S.SLplot)&& S.dt==(tplot-tnow)/365
                xp_mc{tsl,:}=x_mc;
                yp_mc{tsl,:}=y_mc;
                tsl=tsl+1;
                tplot=plot_time(tsl);
            end
            if ~isempty(S.SLplot) && tsl>1
                LName=S.LDBplot(:,2);
                for sl=1:size(xp_mc,1)
                    hp6(sl+mm)=plot(xp_mc{sl,:},yp_mc{sl,:},S.SLplot{sl,3},'linewidth',1.5);
                    LName(end+1)=S.SLplot(sl,2);
                end
                hleg = legend(hp6,LName','Location',S.legendlocation);
                set(hleg,'Box','off','Color','None');
            elseif isempty(S.SLplot) && ~isempty(S.LDBplot)
                hleg = legend(hp6,S.LDBplot(:,2)','Location',S.legendlocation);
                set(hleg,'Box','off','Color','None'); 
            end
             %for shorelines plotting and extracting
           

            %% FORMATTING
            hold off;
            axis equal;
            xlim(S.xlimits);
            ylim(S.ylimits);
            set(gca,'XtickLabel',num2str(get(gca,'Xtick')'/1000,'%2.1f'));
            set(gca,'YtickLabel',num2str(get(gca,'Ytick')'/1000,'%2.1f'));
            xlabel('Easting [km]');
            ylabel('Northing [km]');
            drawnow;
            
            %video
            vi(it+1)=getframe(figure(11));
            pause(0.01);
           
            figstep=round((1/S.fignryear)/S.dt);
            if mod(it,figstep)==0
                %fname=[num2str(it+1000),'.jpg'];
                fname=[num2str(round((it+1)+100000)),'.jpg'];
                setFIGUREproperties(800,600,32);
                if ~exist(S.outputdir)
                    mkdir(S.outputdir);
                end
                
                print('-djpeg','-r300',[S.outputdir,filesep,fname]);
            end
            %         subplot(122)
            %         plot(philoc(xS<1200),yS(xS<1200),'o')
            %         ylim([2000 8000])
            %         title(['section ',num2str(i_mc)])
            %         %xlim([-1e6 1e6])
        end
        S.xyt(it+1).x_mc=x_mc;
        S.xyt(it+1).y_mc=y_mc;

    end
    S.x_mc=x_mc;
    S.y_mc=y_mc;
    S.x_mc0=x_mc0;
    S.y_mc0=y_mc0;
    S.timenum0=timenuma(end);
    % Apply shoreline change due to tide
    if S.tide_interaction
        [S]=update_shoreline(S);
    end
tnow = tnow+ S.dt*365; %dt [year]   %calculate the current time after each time step
end
%print((strcat(S.trform,'_',num2str(S.b/1e4),'_ds0_',num2str(S.ds0),'_phi_',num2str(S.phiw0/pi*180))),'-dpng', '-r300')
%print((strcat(S.trform,'_',num2str(k*100),'_ds0_',num2str(S.ds0),'_phi_',num2str(S.phiw0/pi*180))),'-dpng', '-r300')
%% Plot shorelines at specific dates
            
            if ~isempty(S.SLplot)
                %xp_mc(end+1,:)=x_mc;
                %yp_mc(end+1,:)=y_mc;
                %S.SLplot{end+1,3}='b-'
                %S.SLplot{end+1,2}=''
                L2Name=S.LDBplot(:,2);
                figure(12) 
                for mm=1:size(S.LDBplot,1)
                   
                    LDBplotval = landboundary('read',S.LDBplot{mm,1});
                    hp6(mm)=plot(LDBplotval(:,1)-S.XYoffset(1),LDBplotval(:,2)-S.XYoffset(2),S.LDBplot{mm,3},'linewidth',1.5);
                 hold on
                end
                for sl=1:size(xp_mc,1)
                    hp6(sl+mm)=plot(xp_mc{sl,:},yp_mc{sl,:},S.SLplot{sl,3},'linewidth',1.5);
                    L2Name(end+1)=S.SLplot(sl,2);
                end
                hlegs = legend(hp6,L2Name','Location',S.legendlocation);
                set(hlegs,'Box','off','Color','None');
                plot(x_hard,y_hard,'k','linewidth',2);
                plot(x_mc0,y_mc0,'b','linewidth',2);
                hold off;
                xlim(S.xlimits);
                ylim(S.ylimits);
                set(gca,'XtickLabel',num2str(get(gca,'Xtick')'/1000,'%2.1f'));
                set(gca,'YtickLabel',num2str(get(gca,'Ytick')'/1000,'%2.1f'));
                xlabel('Easting [km]');
                ylabel('Northing [km]');
              % print((strcat(S.trform,'_',num2str(S.b/1e4),'_ds0_',num2str(S.ds0),'_phi_',num2str(S.phiw0/pi*180))),'-dpng', '-r300')
            end
            
            %%videos
            label=strcat(S.trform,'_',S.boundary_condition_start,'_ds0_',num2str(S.ds0),'_tc_',num2str(S.tc));
            FR=25;
            video=VideoWriter(label,'MPEG-4');
            video.FrameRate=FR;
            open(video)
            writeVideo(video,vi)
            close (video)

function [RUNUP]=introduce_runup(RUNUP,TIME,COAST,DUNE,STRUC)
% function [RUNUP]=introduce_runup(RUNUP,TIME,COAST,DUNE,STRUC)
% Interpolate runup conditions SWL, HS, TP and DIR in time
% and along the coast

    if DUNE.used || STRUC.transmission==1
        %% Interpolate conditions in time (all columns for all positions at once)
        x=COAST.x;
        y=COAST.y;

        %% interpolate SWL
        SWL=RUNUP.SWL;
        SWLnow=[];
        for kk=1:length(SWL)
            % interpolate swl at tnow
            eps=1e-6;          
            
            % check how many time instances of the SWL are within the current coastline timestep
            idt=find(SWL(kk).timenum>TIME.tprev & SWL(kk).timenum<=TIME.tnow);
            
            % construct dune timesteps in case a different dt is specified for the dunes than for the coastline
            timenow=TIME.tnow;
            if ~isempty(DUNE.dt) && TIME.it~=0
                dt=min(DUNE.dt*365,TIME.tnow-TIME.tprev);
                timenow=[TIME.tprev+dt:dt:TIME.tnow]';
            end
            
            % interpolate the water-levels
            if length(idt)>1 && min(abs(SWL(kk).timenum-TIME.tnow))<eps && min(abs(SWL(kk).timenum-TIME.tprev))<eps && isempty(DUNE.dt)
                % in case more dense time-series than timesteps of coastline model
                SWLnow(:,kk) = SWL(kk).swl(idt);
            elseif length(SWL(kk).timenum)>1 || TIME.it==0
                % in case similar or larger timestep in data time-series than timesteps of coastline model
                SWLnow(:,kk) = interp1(SWL(kk).timenum,SWL(kk).swl,timenow);
            else
                % in case a single value is specified as the swl
                SWLnow(1,kk) = SWL(kk).swl;
            end
        end
        if RUNUP.nloc==0
            RUNUP.swl=RUNUP.swl(1)*ones(size(COAST.x));
        elseif RUNUP.nloc==1
            % Use uniform swl
            RUNUP.swl=SWLnow.*ones(size(SWLnow,1),length(COAST.x));
        else
            % Interpolate conditions along the coast
            method='weighted_distance';
            xw=RUNUP.x;
            yw=RUNUP.y;
            RUNUP.swl=get_interpolation_on_grid(method,x,y,xw,yw,SWLnow);
        end
    end
    if DUNE.used
        %% interpolate WVD
        WVD=RUNUP.WVD;
        for kk=1:length(WVD)
            % interpolate swl at tnow
            eps=1e-6;          
            
            % check how many time instances of the SWL are within the current coastline timestep
            idt=find(WVD(kk).timenum>TIME.tprev & WVD(kk).timenum<=TIME.tnow);
            
            % construct dune timesteps in case a different dt is specified for the dunes than for the coastline
            timenow=TIME.tnow;
            if ~isempty(DUNE.dt) && TIME.it~=0
                dt=min(DUNE.dt*365,TIME.tnow-TIME.tprev);
                timenow=[TIME.tprev+dt:dt:TIME.tnow]';
            end
            
            % interpolate the wave conditions
            if length(idt)>1 && min(abs(WVD(kk).timenum-TIME.tnow))<eps && min(abs(WVD(kk).timenum-TIME.tprev))<eps && isempty(DUNE.dt)
                % in case more dense time-series than timesteps of coastline model
                HSnow(:,kk) = WVD(kk).Hs(idt);
                TPnow(:,kk) = WVD(kk).Tp(idt);
                DIRnow(:,kk)= WVD(kk).Dir(idt);
            elseif length(WVD(kk).timenum)>1 || TIME.it==0
                % in case similar or larger timestep in data time-series than timesteps of coastline model
                HSnow(:,kk) = interp1(WVD(kk).timenum,WVD(kk).Hs,timenow);
                TPnow(:,kk) = interp1(WVD(kk).timenum,WVD(kk).Tp,timenow);
                sinDIRnow   = interp1(WVD(kk).timenum,sind(WVD(kk).Dir),timenow);
                cosDIRnow   = interp1(WVD(kk).timenum,cosd(WVD(kk).Dir),timenow);
                DIRnow(:,kk)= mod(atan2d(sinDIRnow,cosDIRnow),360);
            else
                % in case a single value is specified as the WVD
                HSnow(1,kk) = WVD(kk).Hs;
                TPnow(1,kk) = WVD(kk).Tp;
                DIRnow(1,kk)= WVD(kk).Dir;
            end
        end
        if RUNUP.nlocw==0
            RUNUP.Hs=RUNUP.Hs(1)*ones(size(COAST.x));
            RUNUP.Tp=RUNUP.Tp(1)*ones(size(COAST.x));
            RUNUP.Dir=RUNUP.Dir(1)*ones(size(COAST.x));
        elseif RUNUP.nlocw==1
            RUNUP.Hs=HSnow*ones(size(COAST.x));
            RUNUP.Tp=TPnow*ones(size(COAST.x));
            RUNUP.Dir=DIRnow*ones(size(COAST.x));
        else
            % Interpolate conditions along the coast
            method='weighted_distance';
            x=COAST.x;
            y=COAST.y;
            xw=RUNUP.xw;
            yw=RUNUP.yw;
            var1=[];
            var2=[];
            var1.Hs=HSnow;
            var1.Tp=TPnow;
            var2.Dir=DIRnow;
            [var1i,var2i]=get_interpolation_on_grid(method,x,y,xw,yw,var1,var2);
            RUNUP.Hs=var1i.Hs;
            RUNUP.Tp=var1i.Tp;
            RUNUP.Dir=var2i.Dir;
        end
    end
end
function [S,O]=ShorelineS(S0)
% MODEL : ShorelineS
% 
% This model computes shoreline changes as a result of gradients in alongshore
% sediment transport for arbitrary shaped coastlines.
% 
% INPUT: 
%     S       data structure wiht fields:
%              .<properties>
% 
% created by:   J.A. Roelvink (2016-present) - IHE Delft
% extended by:  B.J.A. Huisman (2017-present) - Deltares
% 
% additions to the code:
%  A.M. Elghandour (2018-present) - IHE Delft
%  J. Reyns (2018-present) - IHE Delft, Deltares
%  M.E. Ghonim (2019) - IHE Delft
%  C.M. Mudde (2019) - TU-Delft
%  B. Perry (2022) - Flinders University
%  A. de Beer (2023) - Deltares 
%  A. de Bakker (2023) - Deltares
%  K. Trouw (2023) - Flanders Hydraulics
% 
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2020 IHE Delft & Deltares
%
%       Dano Roelvink
%       d.roelvink@un-ihe.org
%       Westvest 7
%       2611AX Delft
%
%       Bas Huisman
%       bas.huisman@deltares.nl
%       Boussinesqweg 1
%       2629HV Delft
%
%   This library is free software: you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation, either
%   version 2.1 of the License, or (at your option) any later version.
%
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses>
%   --------------------------------------------------------------------

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DEFAULT INPUT PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [S]       = initialize_defaultvalues(S0);
    [S]       = initialize_randomgenerator(S);
    [TIME]    = initialize_time(S);
    [COAST]   = prepare_coastline(S);
    [DUNE]    = prepare_dunes(S,COAST);
    [CC]      = prepare_climatechange(S,TIME);
    [WAVE]    = prepare_waveconditions(S,TIME);
    [RUNUP]   = prepare_runupconditions(S,TIME);
    [WIND]    = prepare_windconditions(S,TIME);   
    [TIDE]    = prepare_tide(S);
    [STRUC]   = prepare_structures(S,COAST);
    [NOUR]    = prepare_nourishment(S,COAST,STRUC);
    [FNOUR]   = prepare_fnourishment(S);
    [TRANSP]  = prepare_transport(S);
    [MUD]     = prepare_mudcoast(S);
    [SPIT]    = prepare_spit(S);
    [CHANNEL] = prepare_channel(S);
    [DELTA]   = prepare_delta(S);
    [BATHY]   = initialize_bathyupdate(S);
    [FORMAT]  = initialize_plot(S,COAST);    
    [O,P,V]   = initialize_output(S,TIME,COAST,DUNE,MUD);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Loop over time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('  Loop over time \n'); 
    while TIME.tnow<TIME.tend || (TIME.tnow==TIME.tend && TIME.tc==0)
        TIME.it=TIME.it+1; 
        TIME.nt=TIME.it; 
        
        if S.debug==2 
            if ~isoctave, warning off, end 
            save('debug.mat'); 
            if ~isoctave, warning on, end
        end 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PHASE 0 : GRID                                        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [COAST,STRUC,GROYNE]=prepare_grid_groyne(COAST,STRUC,S.yesplot);   
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PHASE 1 : TRANSPORT                                        %%
        %% loop over coastline sections                               %%
        %% compute sediment transport                                 %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i_mc=1:COAST.n_mc 
            
            %% Alongshore coordinate s; regridding where necessary
            [COAST]=make_sgrid_mc(COAST,TIME,i_mc);
            
            %% Interpolate or update all properties on coast points
            [COAST]=interpolate_props(COAST,DUNE,MUD,TIME);
            
            %% Determine climate change effects
            [CC]=introduce_climatechange(CC,TIME,i_mc);
            
            %% Interpolate wave conditions along the coast
            [WAVE]=introduce_wave(WAVE,TIME,COAST,CC);
                        
            %% Interpolate runup conditions along the coast
            [RUNUP]=introduce_runup(RUNUP,TIME,COAST,DUNE,STRUC,WAVE);
            
            %% Interpolate tide components along the coast
            [TIDE]=interpolate_tide(TIDE,COAST);
            
            %% Get coastline orientation PHIc
            [COAST]=get_coastline_orientation(COAST);
            
            %% Get foreshore orientation PHIf
            [COAST]=get_foreshore_orientation(COAST);
            
            %% Get active profile 
            [COAST]=get_active_profile(COAST);
            
            %% Get refracted waves
            % From here on, WAVE.PHItdp and WAVE.HStdp are always given at the toe of the
            % dynamic profile (TDP) and are always with a size of 1 by n
            [WAVE]=wave_refraction(WAVE,COAST,TRANSP);
            
            %% Wave diffraction
            [STRUC,WAVE]=wave_diffraction(STRUC,COAST,WAVE,GROYNE,RUNUP);
            
            %% Introduce the wind conditions to the domain
            [WIND]=introduce_wind(WIND,TIME,WAVE,COAST);
            
            %% Relative angle of the offshore and nearshore waves 
            WAVE.dPHIo=atan2d(sind(COAST.PHIc-WAVE.PHIo),cosd(COAST.PHIc-WAVE.PHIo));
            WAVE.dPHItdp=atan2d(sind(COAST.PHIcs-WAVE.PHItdp),cosd(COAST.PHIcs-WAVE.PHItdp));   
            
            %% Wave height in the nearshore (due to refraction and shoaling)
            [WAVE]=wave_breakingheight(WAVE,TRANSP);
            
            %% Nearshore wave direction at point of breaking
            WAVE.PHIbr=mod(COAST.PHIc-WAVE.dPHIbr,360); 
            
            %% Introduce permeable structures
            [WAVE]=introduce_perm_structures(STRUC,WAVE,COAST);
            
            %% Critical angles for transport -> WAVE.dPHIcrit & TRANSP.QSmax
            [WAVE,TRANSP]=wave_angles(COAST,WAVE,TIDE,TRANSP,STRUC);
            
            %% Longshore Transport
            [TRANSP]=transport(TRANSP,WAVE,TIDE,STRUC);
            
            %% River discharges for mud transport
            [MUD]=get_riverdischarges(TIME,COAST,TRANSP,MUD);
            
            %% Mud transport
            [TRANSP,COAST]=transport_mud(COAST,TRANSP,WAVE,WIND,MUD);
            
            %% Shadowing effect on Transport
            [TRANSP,WAVE]=transport_shadow_treat(COAST,STRUC,WAVE,TRANSP);
            
            %% Revetments
            [TRANSP,COAST]=transport_revetment(COAST,STRUC,TRANSP,TIME,WAVE);  
            
            %% Smoothen shoreline angles
            [COAST,TRANSP]=get_smoothangles(COAST,TRANSP);
            
            %% Upwind correction for high-angle -> using wave angle at the nearshore location ('tdp')
            [TRANSP]=get_upwindcorrection(COAST,WAVE,TRANSP);          
            
            %% Sand Bypassing and Transmission (GROYNE.QS)
            [GROYNE]=transport_bypass(TRANSP,WAVE,COAST,STRUC,GROYNE,TIDE);
            
            %% Cross-shore flux from/to dune
            [DUNE,COAST]=dune_flux(COAST,DUNE,WIND,RUNUP,TRANSP,TIME,CC);
            
            %% Boundary condition
            [TRANSP]=transport_boundary_condition(TRANSP,COAST,GROYNE);      
            
            %% Adaptive time step based on transport of separate coastline sections
            [TIME]=get_timestep(TIME,COAST,TRANSP,WAVE);
            
            %% Collect QS, s and WAVE.PHItdp in QS_mc, WAVE.PHIo_mc and s_mc -> stores the data of this coastline section in 'mc'
            [COAST,WAVE,TRANSP]=collect_variables(COAST,WAVE,TRANSP,DUNE,MUD,GROYNE,TIME,i_mc);
            
            %% debugging plot for QS and SPHI
            plot_debug(S.debug,COAST,WAVE,TRANSP); % only used when S.debug=1;
            
        end
        TIME.adt_record(TIME.it+1)=TIME.adt;
        
        %% Final timestep for coastline update 
        [TIME,WAVE]=get_coastlineupdate_timestep(TIME, BATHY, WAVE, FORMAT);    
        
        if S.debug==2
            if ~isoctave, warning off, end
            save('debug2.mat');
            if ~isoctave, warning on, end
        end    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PHASE 2 : COASTLINE CHANGE                                 %%
        %% Compute coastline and dune change                          %%
        %% Add nourishments to the coast                              %%
        %% Re-connect coastal sections at groynes                     %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i_mc=1:COAST.n_mc
            
            %% Retrieve data of individual segment (i_mc)
            [COAST,WAVE,TRANSP]=get_segmentdata(COAST,WAVE,TRANSP,DUNE,MUD,i_mc);
            
            %% Nourishment
            [NOUR]=get_nourishments(TIME,COAST,NOUR);
            
             %% Shoreface nourishment 
            [FNOUR]=get_fnourishment(TIME,COAST,FNOUR,WAVE,i_mc);
            
            %% Coastline change 
            [COAST,GROYNE,STRUC,TRANSP]=coastline_change(COAST,WAVE,TRANSP,DUNE,MUD,STRUC,GROYNE,TIME,NOUR,FNOUR,CC);
            
        end % loop over sections
        
        %% Update the dune and mud properties
        [COAST,DUNE,MUD] = update_props(COAST,DUNE,MUD);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PHASE 3 : OTHER PROCESSES                                  %%
        %% Overwash process                                           %%
        %% Reconnect the section separated by groines                 %%
        %% Move channels                                              %%
        %% Splitting and Merging coastlines                           %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Overwash process
        [COAST,SPIT]=find_overwash_mc(COAST,WAVE,SPIT,STRUC,TRANSP,DUNE,MUD,TIME);
        
        %% Move channel
        [CHANNEL,COAST]=move_channel(CHANNEL,COAST,TIME);
        
        %% reconnect coastlines through groynes
        [COAST]=get_reconnectedgroynes(COAST,GROYNE); 
        
        %% merge multiple coastline sections -> only x_mc and y_mc are currently merged.
        for i_mc=1:COAST.n_mc
            [COAST]=merge_coastlines(COAST,i_mc);  
        end
        [COAST]=merge_coastlines_mc(COAST);

        %% make transport points xq_mc and yq_mc
        [COAST]=get_transportpoints(COAST,1:COAST.n_mc);
        
        %% clean up redundant NaNs
        [COAST]=cleanup_nans(COAST);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PHASE 4 : PLOTTING AND STORING COASTLINES                  %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        if O.netcdf==1 % write to 'shorelines_output.nc'
            [O,P]=save_shorelines_netcdf(O,P,TIME,COAST,WAVE,TRANSP,DUNE,MUD);
        else % write to 'output.mat'
            [O,P]=save_shorelines(O,P,TIME,COAST,WAVE,TRANSP,STRUC,NOUR,FNOUR,GROYNE,DUNE,MUD); 
        end

        %% Plotting                                     
        [V,FORMAT,TIME]=plot_coast(CHANNEL,STRUC,COAST,DUNE,WAVE,TIME,TRANSP,FORMAT,V,FNOUR);
        [FORMAT]=plot_profiles(FORMAT,TIME,DUNE,O);
        
        %% Apply shoreline change due to tide           
        [FORMAT,COAST]=update_shoreline(S,COAST,FORMAT);
        
        %% Bathymetry update & Plot                              
        [S,BATHY]=update_bathy(S,BATHY,TIME,COAST);
        
        %% Next time step + Status update on screen
        [TIME]=get_nexttimestep(TIME,WAVE,DUNE,WIND);
    end
       
    %% Extract shorelines coordinates and plotting figures at specific dates (x_mc0 and y_mc0)
    fprintf('  Post-process results \n');
    extract_shoreline(S,STRUC,COAST,FORMAT);
    
    %% Capture video frames
    if S.video==1
        make_video(S,V);
    end
end
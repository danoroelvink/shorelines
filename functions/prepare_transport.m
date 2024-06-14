function [TRANSP]=prepare_transport(S)
% function [TRANSP]=prepare_transport(S)
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
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses
%   --------------------------------------------------------------------

    fprintf('  Prepare transport \n');
    TRANSP=struct;
    TRANSP.trform=S.trform;                                                            % switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
    TRANSP.b=S.b;                                                                      % CERC : coeff in simple cerc formula
    TRANSP.qscal=S.qscal;                                                              % Calibration factor of the transport (works for all transport formulas)
    TRANSP.d50=S.d50;                                                                  % KAMP & MILH & VR14 : median grain diameter [m]
    TRANSP.d90=S.d90;                                                                  % KAMP & MILH & VR14 : d90 grain diameter [m]
    TRANSP.porosity=S.porosity;                                                        % KAMP & MILH & VR14 : TRANSP.porosity (typically 0.4) [-]
    TRANSP.tanbeta=S.tanbeta;                                                          % KAMP & MILH & VR14 : mean bed slope [ratio 1/slope]
    TRANSP.ks=S.ks;
    TRANSP.rhos=S.rhos;                                                                % KAMP & MILH & VR14 : density of sand [kg/m3]
    TRANSP.rhow=S.rhow;                                                                % KAMP & MILH & VR14 : density of water [kg/m3]
    TRANSP.g=S.g;                                                                      % KAMP & MILH & VR14 : gravitational acceleration [m2/s]
    TRANSP.Cf=S.Cf;
    TRANSP.alpha=S.alpha;                                                              % KAMP & MILH & VR14 : calibration factor for point of breaking (TRANSP.alpha = 1.8 for Egmond data)
    TRANSP.gamma=S.gamma;                                                              % KAMP & MILH & VR14 : breaking coefficient (Hs/h) with 5% breaking waves
    TRANSP.Pswell=S.Pswell;                                                            % VR14 : Percentage swell (between 0 - 100) [-]
    TRANSP.crit=S.crit;                                                                % stability criterion (not active)
    TRANSP.Aw=S.Aw;                                                                    % factor for determining depth of closure at bypassing groyne (1.27 if time series is used) This value is used by default.
    TRANSP.bypasscontractionfactor=S.bypasscontractionfactor;                          % the maximum transport when the coastline is at the end of the structure (always >=1). Setting the bypass fraction larger than 1 means that the accretion does not go to the tip of the structure. 
    TRANSP.relaxationlength=S.relaxationlength;                                        % length over which transport decelerates in meters, which adds inertia to the longshore current. It scales linearly with the wave height below 1m waves.
    if isempty(S.WVCfile)
        TRANSP.Aw=S.Awfixedhs;                                                         % factor for determining depth of closure at bypassing groyne ir a representative Hs is used instead of a climate or timeseries. This value is used instead of 'Aw' if S.WVCfile is empty. 
    end
    TRANSP.twopoints=S.twopoints;
    TRANSP.crit_width=S.crit_width;
    TRANSP.suppress_highangle=S.suppress_highangle;                                    % switch 0/1 to disable the high-angle instabilities by limiting the transport angle to the critical high-angle orientation (when it is set at 1)
    % boundary conditions
    % {'Closed',e.g. 0 or 9000 m3/yr}, {'Neumann',dummy},{'Fixed',dummy},{'Angleconstant',empty to use at t0 or specified value e.g. 321�N}; 
    if ischar(S.boundary_condition_start)
        TRANSP.boundary_condition_start={S.boundary_condition_start,nan};                     % boundary condition 'Closed', 'Neumann','Fixed','Angleconstant'
    else
        TRANSP.boundary_condition_start=S.boundary_condition_start;
    end
    if ischar(S.boundary_condition_end) 
        TRANSP.boundary_condition_end={S.boundary_condition_end,nan};                       % boundary condition 'Closed', 'Neumann','Fixed','Angleconstant'
    else
        TRANSP.boundary_condition_end=S.boundary_condition_end;
    end 
    
    
    %% Sediment limitation ('virtual revetment')
    TRANSP.x_sedlim=[];     % x-coordinates for designated sediment limitation areas along the coast
    TRANSP.y_sedlim=[];     % y-coordinates for designated sediment limitation areas along the coast
    %TRANSP.xc_sedlim=[];    % cross-shore distance of sediment limiter line w.r.t. initial coastline
    TRANSP.width_sedlim=[];    % cross-shore distance with 100% transport of sediment w.r.t. position of sediment limiter line
    TRANSP.sedlim=S.sedlim;
    if TRANSP.sedlim
        %  add sediment limitations interactively with the graphical user interface
        if strcmpi(S.LDBsedlim,'manual') || strcmpi(S.LDBsedlim,'interactive') 
            figure(11);
            axis equal;
            xl=xlim;yl=ylim;
            htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add sediment limitation element (LMB); Next sediment limitation element (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
            [x_sedlim,y_sedlim]=select_multi_polygon('k');
            set(htxt2,'Visible','off');
            TRANSP.x_sedlim=x_sedlim;
            TRANSP.y_sedlim=y_sedlim;
            TRANSP.width_sedlim=repmat(TRANSP.width_sedlim,size(TRANSP.x_sedlim));
            
        %  try to load a file with sediment limitation 
        elseif ~isempty(S.LDBsedlim)
            xy_sedlim=load(S.LDBsedlim);
            TRANSP.x_sedlim=xy_sedlim(:,1)'-S.XYoffset(1);
            TRANSP.y_sedlim=xy_sedlim(:,2)'-S.XYoffset(2);
            TRANSP.width_sedlim=xy_sedlim(:,3)';
            
        % add a sediment limitation by directly specifying the x_sedlim and y_sedlim
        elseif ~isempty(S.x_sedlim) && ~isempty(S.y_sedlim)
            TRANSP.x_sedlim=S.x_sedlim(:)'-S.XYoffset(1);
            TRANSP.y_sedlim=S.y_sedlim(:)'-S.XYoffset(2);
            TRANSP.width_sedlim=S.width_sedlim(:)';
            if length(TRANSP.width_sedlim)<length(TRANSP.x_sedlim)
                width_sedlim=TRANSP.width_sedlim(~isnan(TRANSP.width_sedlim));
                % in case TRANSP.width_sedlim is defined as [value_for_section_1, value_for_section_2, ...]
                idnan=find(isnan(TRANSP.x_sedlim));
                if length(width_sedlim)==length(idnan)+1
                    ids=[0,idnan(:)',length(TRANSP.x_sedlim)+1];
                    TRANSP.width_sedlim=nan(size(TRANSP.x_sedlim));
                    jj=1;
                    for ii=1:length(TRANSP.width_sedlim)
                        if ~isnan(TRANSP.x_sedlim(ii))
                            TRANSP.width_sedlim(ii)=width_sedlim(jj);
                        else
                            jj=jj+1;
                        end
                    end
                else
                    % in case TRANSP.width_sedlim is defined as a single value for all sections
                    TRANSP.width_sedlim=repmat(TRANSP.width_sedlim(1),size(TRANSP.x_sedlim));
                end
            end
            
        % no sediment limitations
        else
            TRANSP.sedlim=0;
            TRANSP.x_sedlim=[];
            TRANSP.y_sedlim=[];
            TRANSP.width_sedlim=[];
        end
    end
    
    % tidy up the sediment limitation without nans at the start and end
    idnotnan=find(~isnan(TRANSP.x_sedlim));
    if ~isempty(TRANSP.x_sedlim)
        TRANSP.x_sedlim=TRANSP.x_sedlim(idnotnan(1):idnotnan(end));
        TRANSP.y_sedlim=TRANSP.y_sedlim(idnotnan(1):idnotnan(end));
        TRANSP.width_sedlim=TRANSP.width_sedlim(idnotnan(1):idnotnan(end));
        idnan=find(isnan(TRANSP.x_sedlim));
        iduse=setdiff([1:length(TRANSP.x_sedlim)],idnan(diff(idnan)==1));
        TRANSP.x_sedlim=TRANSP.x_sedlim(iduse);
        TRANSP.y_sedlim=TRANSP.y_sedlim(iduse);
        TRANSP.width_sedlim=TRANSP.width_sedlim(iduse);
    end
    
    %% write logfile
    % struct2log(TRANSP,'TRANSP','a');
       
end

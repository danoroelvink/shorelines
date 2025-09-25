function [TRANSP]=prepare_transport(S)
% function [TRANSP]=prepare_transport(S)
%
% The sediment transport is initialized in this function. 
% A data-structure TRANSP is created, which is used throughout the computation. 
% 
% INPUT:
%   S
%     .trform                   : switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
%     .b                        : CERC : coeff in simple cerc formula
%     .qscal                    : Calibration factor of the transport (works for all transport formulas)
%     .d50                      : median grain diameter [m]
%     .d90                      : d90 grain diameter [m]
%     .porosity                 : TRANSP.porosity (typically 0.4) [-]
%     .tanbeta                  : mean bed slope [ratio 1/slope]
%     .tanbetasetup             : default value for bed slope effect of the water-level setup driven currents impact on sediment transport (dHs/ds), with 1 as default for a very small impact on transport (for a steep slope)
%     .ks                       : roughness parameter
%     .rhos                     : density of sand [kg/m3]
%     .rhow                     : density of water [kg/m3]
%     .g                        : gravitational acceleration [m2/s]
%     .cf                       : roughness factor
%     .alpha                    : calibration factor for point of breaking (TRANSP.alpha = 1.8 for Egmond data)
%     .gamma                    : breaking coefficient (Hs/h) with 5% breaking waves
%     .pswell                   : VR14 : Percentage swell (between 0 - 100) [-]
%     .aw                       : factor for determining depth of closure at bypassing groyne (1.27 if time series is used) This value is used by default.
%     .bypasscontrfac           : the maximum transport when the coastline is at the end of the structure (always >=1). Setting the bypass fraction larger than 1 means that the accretion does not go to the tip of the structure. 
%     .relaxationlength         : length over which transport decelerates in meters, which adds inertia to the longshore current. It scales linearly with the wave height below 1m waves.
%     .aw                       : factor for determining depth of closure at bypassing groyne ir a representative Hs is used instead of a climate or timeseries. This value is used instead of 'aw' if S.wvcfile is empty. 
%     .twopoints                : approach for dealing with high-angle instabilities
%     .critwidth                : critical width of barriers for overwash
%     .suppresshighangle        : switch 0/1 to disable the high-angle instabilities by limiting the transport angle to the critical high-angle orientation (when it is set at 1)
%     .boundaryconditionstart   : left boundary condition of the model (e.g. 'Fixed')
%     .boundaryconditionend     : right boundary condition of the model (e.g. {'Angleconstant',217} )
%     .sedlim                   : switch for sediment limitation 
%     .ldbsedlim   (option 1)   : filename for region with limited sediment availability
%     .xsedlim     (option 2)   : x-coordinates of locations with limited sediment [m]
%     .ysedlim     (option 2)   : y-coordinates of locations with limited sediment [m]
%     .widthsedlim              : cross-shore distance with 100% transport of sediment w.r.t. position of sediment limiter line [m]
%     .xyoffset                 : x and y offset used for plotting [1x2] in [m]
% 
% OUTPUT:
%   TRANSP
%     .trform                   : switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
%     .b                        : CERC : coeff in simple cerc formula
%     .qscal                    : Calibration factor of the transport (works for all transport formulas)
%     .d50                      : median grain diameter [m]
%     .d90                      : d90 grain diameter [m]
%     .porosity                 : TRANSP.porosity (typically 0.4) [-]
%     .tanbeta                  : mean bed slope [ratio 1/slope]
%     .tanbetasetup             : default value for bed slope effect of the water-level setup driven currents impact on sediment transport (dHs/ds), with 1 as default for a very small impact on transport (for a steep slope)
%     .ks                       : roughness parameter
%     .rhos                     : density of sand [kg/m3]
%     .rhow                     : density of water [kg/m3]
%     .g                        : gravitational acceleration [m2/s]
%     .cf                       : roughness factor
%     .alpha                    : calibration factor for point of breaking (TRANSP.alpha = 1.8 for Egmond data)
%     .gamma                    : breaking coefficient (Hs/h) with 5% breaking waves
%     .pswell                   : VR14 : Percentage swell (between 0 - 100) [-]
%     .aw                       : factor for determining depth of closure at bypassing groyne (1.27 if time series is used) This value is used by default.
%     .bypasscontrfac           : the maximum transport when the coastline is at the end of the structure (always >=1). Setting the bypass fraction larger than 1 means that the accretion does not go to the tip of the structure. 
%     .relaxationlength         : length over which transport decelerates in meters, which adds inertia to the longshore current. It scales linearly with the wave height below 1m waves.
%     .aw                       : factor for determining depth of closure at bypassing groyne ir a representative Hs is used instead of a climate or timeseries. This value is used instead of 'aw' if S.wvcfile is empty. 
%     .twopoints                : approach for dealing with high-angle instabilities
%     .critwidth                : critical width of barriers for overwash
%     .suppresshighangle        : switch 0/1 to disable the high-angle instabilities by limiting the transport angle to the critical high-angle orientation (when it is set at 1)
%     .boundaryconditionstart   : left boundary condition of the model (e.g. 'Fixed')
%     .boundaryconditionend     : right boundary condition of the model (e.g. {'Angleconstant',217} )
%     .sedlim                   : switch for sediment limitation 
%     .xsedlim                  : x-coordinates of locations with limited sediment [m]
%     .ysedlim                  : y-coordinates of locations with limited sediment [m]
%     .nsedlim                  : number of elements with limited sediment
%     .widthsedlim              : cross-shore distance with 100% transport of sediment w.r.t. position of sediment limiter line [m]
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

    fprintf('  Prepare transport \n');
    TRANSP=struct;
    TRANSP.trform=S.trform;                                                            % switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
    TRANSP.b=S.b;                                                                      % CERC : coeff in simple cerc formula
    TRANSP.qscal=S.qscal;                                                              % Calibration factor of the transport (works for all transport formulas)
    TRANSP.d50=S.d50;                                                                  % median grain diameter [m]
    TRANSP.d90=S.d90;                                                                  % d90 grain diameter [m]
    TRANSP.porosity=S.porosity;                                                        % TRANSP.porosity (typically 0.4) [-]
    TRANSP.tanbeta=S.tanbeta;                                                          % mean bed slope used for computing the transport in formulas [ratio 1/slope]
    TRANSP.tanbetasetup=S.tanbetasetup;                                                % bed slope value, scaling the effect of the water-level setup driven currents impact on sediment transport (dHs/ds), with 1 as default for a very small impact on transport (for a steep slope)
    TRANSP.ks=S.ks;
    TRANSP.rhos=S.rhos;                                                                % density of sand [kg/m3]
    TRANSP.rhow=S.rhow;                                                                % density of water [kg/m3]
    TRANSP.g=S.g;                                                                      % gravitational acceleration [m2/s]
    TRANSP.cf=S.cf;                                                                    % roughness factor [m]
    TRANSP.n=S.n;                                                                      % manning coefficient for wave induced currents in the tide module [s/m^(1/3)]
    TRANSP.alpha=S.alpha;                                                              % calibration factor for point of breaking (TRANSP.alpha = 1.8 for Egmond data)
    TRANSP.gamma=S.gamma;                                                              % breaking coefficient (Hs/h) with 5% breaking waves
    TRANSP.pswell=S.pswell;                                                            % VR14 : Percentage swell (between 0 - 100) [-]
    TRANSP.aw=S.aw;                                                                    % factor for determining depth of closure at bypassing groyne (1.27 if time series is used) This value is used by default.
    if isempty(S.wvcfile)
        TRANSP.aw=S.awfixedhs;                                                         % factor for determining depth of closure at bypassing groyne ir a representative Hs is used instead of a climate or timeseries. This value is used instead of 'aw' if S.wvcfile is empty. 
    end
    TRANSP.relaxationlength=S.relaxationlength;                                        % length over which transport decelerates in meters, which adds inertia to the longshore current. It scales linearly with the wave height below 1m waves.
    TRANSP.twopoints=S.twopoints;
    TRANSP.critwidth=S.critwidth;
    TRANSP.suppresshighangle=S.suppresshighangle;                                      % switch 0/1 to disable the high-angle instabilities by limiting the transport angle to the critical high-angle orientation (when it is set at 1)
    TRANSP.acal=S.acal;                                                                % calibration factor for the Sulsby Van Rijn formulation in the tide module
    TRANSP.submerged=S.submerged;
    
    % boundary conditions
    % {'Closed',e.g. 0 or 9000 m3/yr}, {'Neumann',dummy},{'Fixed',dummy},{'Angleconstant',empty to use at t0 or specified value e.g. 321�N}; 
    if ischar(S.boundaryconditionstart)
        TRANSP.boundaryconditionstart={S.boundaryconditionstart,nan};                     % boundary condition 'Closed', 'Neumann','Fixed','Angleconstant'
    else
        TRANSP.boundaryconditionstart=S.boundaryconditionstart;
    end
    if ischar(S.boundaryconditionend) 
        TRANSP.boundaryconditionend={S.boundaryconditionend,nan};                       % boundary condition 'Closed', 'Neumann','Fixed','Angleconstant'
    else
        TRANSP.boundaryconditionend=S.boundaryconditionend;
    end 
    
    %% Sediment limitation ('virtual revetment')
    TRANSP.xsedlim=[];     % x-coordinates for designated sediment limitation areas along the coast
    TRANSP.ysedlim=[];     % y-coordinates for designated sediment limitation areas along the coast
    %TRANSP.xc_sedlim=[];    % cross-shore distance of sediment limiter line w.r.t. initial coastline
    TRANSP.widthsedlim=[];    % cross-shore distance with 100% transport of sediment w.r.t. position of sediment limiter line
    TRANSP.sedlim=S.sedlim;
    if TRANSP.sedlim
        %  add sediment limitations interactively with the graphical user interface
        if strcmpi(S.ldbsedlim,'manual') || strcmpi(S.ldbsedlim,'interactive') 
            figure(11);
            axis equal;
            xl=xlim;yl=ylim;
            htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add sediment limitation element (LMB); Next sediment limitation element (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
            [xsedlim,ysedlim]=select_multi_polygon('k');
            set(htxt2,'Visible','off');
            TRANSP.xsedlim=xsedlim;
            TRANSP.ysedlim=ysedlim;
            TRANSP.widthsedlim=repmat(TRANSP.widthsedlim,size(TRANSP.xsedlim));
            
        %  try to load a file with sediment limitation 
        elseif ~isempty(S.ldbsedlim)
            xysedlim=load(S.ldbsedlim);
            TRANSP.xsedlim=xysedlim(:,1)'-S.xyoffset(1);
            TRANSP.ysedlim=xysedlim(:,2)'-S.xyoffset(2);
            TRANSP.widthsedlim=xysedlim(:,3)';
            
        % add a sediment limitation by directly specifying the xsedlim and ysedlim
        elseif ~isempty(S.xsedlim) && ~isempty(S.ysedlim)
            TRANSP.xsedlim=S.xsedlim(:)'-S.xyoffset(1);
            TRANSP.ysedlim=S.ysedlim(:)'-S.xyoffset(2);
            TRANSP.widthsedlim=S.widthsedlim(:)';
            if length(TRANSP.widthsedlim)<length(TRANSP.xsedlim)
                widthsedlim=TRANSP.widthsedlim(~isnan(TRANSP.widthsedlim));
                % in case TRANSP.widthsedlim is defined as [value_for_section_1, value_for_section_2, ...]
                idnan=find(isnan(TRANSP.xsedlim));
                if length(widthsedlim)==length(idnan)+1
                    ids=[0,idnan(:)',length(TRANSP.xsedlim)+1];
                    TRANSP.widthsedlim=nan(size(TRANSP.xsedlim));
                    jj=1;
                    for ii=1:length(TRANSP.widthsedlim)
                        if ~isnan(TRANSP.xsedlim(ii))
                            TRANSP.widthsedlim(ii)=widthsedlim(jj);
                        else
                            jj=jj+1;
                        end
                    end
                else
                    % in case TRANSP.widthsedlim is defined as a single value for all sections
                    TRANSP.widthsedlim=repmat(TRANSP.widthsedlim(1),size(TRANSP.xsedlim));
                end
            end
            
        % no sediment limitations
        else
            TRANSP.sedlim=0;
            TRANSP.xsedlim=[];
            TRANSP.ysedlim=[];
            TRANSP.widthsedlim=[];
        end
        TRANSP.nsedlim=min(sum(isnan(TRANSP.xsedlim))+1,length(TRANSP.xsedlim));
    end
    
    % tidy up the sediment limitation without nans at the start and end
    idnotnan=find(~isnan(TRANSP.xsedlim));
    if ~isempty(TRANSP.xsedlim)
        TRANSP.xsedlim=TRANSP.xsedlim(idnotnan(1):idnotnan(end));
        TRANSP.ysedlim=TRANSP.ysedlim(idnotnan(1):idnotnan(end));
        TRANSP.widthsedlim=TRANSP.widthsedlim(idnotnan(1):idnotnan(end));
        idnan=find(isnan(TRANSP.xsedlim));
        iduse=setdiff([1:length(TRANSP.xsedlim)],idnan(diff(idnan)==1));
        TRANSP.xsedlim=TRANSP.xsedlim(iduse);
        TRANSP.ysedlim=TRANSP.ysedlim(iduse);
        TRANSP.widthsedlim=TRANSP.widthsedlim(iduse);
    end
    
    %% write logfile
    % struct2log(TRANSP,'TRANSP','a');
       
end

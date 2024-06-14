function make_tidefile(xystat,cZWL,vmean,xsdst,dps,varargin)
% FUNCTION make_tidefile(xystat,wl) creates the tidefile for ShorelineS
%   This function computes M2 and M4 tidal components, given any 
%   timeseries of waterlevels. 
%   Then it creates the tidefile (in the present working directory) as 
%   expected by the tide-module of ShorelineS
%
% INPUT:
%   xystat: [N x 2] matrix containing x- and y-coordinates of N tide 
%           stations (measured or modelled). Minimal 2 stations should be 
%           defined, and ordered in such a way that the sea is on the 
%           lefthand side of the line (conform ShorelineS coastline
%           definition).
%   cZWL:   {N x 1}[nT x 2] cell-array containing N waterlevel timeseries,
%           with timestamp and value in the columns of the matrix. The
%           number of timestamps (nT) can be different for each station (N)
%   vmean:  [N x 1] vector containing the tidally averaged (residual) 
%           alongshore velocity in the stations (if only a scalar is given,
%           it is repeated N times)
%   xsdst:  [N x 1] cross-shore distance (m)
%   dps:    [N x 1] bed level (m MSL)
%   <varargin> keyword-value pairs
%   latitude: switch for nodal correction, activated by vector "lat":
%             [N x 1] vector containing the latitude of the stations 
%             (if only a scalar is given, it is repeated N times)
%   ssSwitch: switch for the computation of the alongshore surface slope:
%             either "harmonic" or "1D_analytical" (default)
%   Cf:       bed friction coefficient (-)
%   hmin:     minimum water depth (m)
%   tidefile: [char-array] name of outputfile
%   
% NOTE: The number of timestamps (nT) can be different for each station 
%       (N): e.g. different time intervals, but it is advisable that they
%       span the same period. Time should be expressed as MATLAB datenum.
%       Having NaN's in the timeseries should be avoided.
%
% OUTPUT:
%   tidefile containing xstat, ystat, etaM2, etaM4, detadsM2, detadsM4,
%   phiM2, phiM4, kM2, kM4, ss for every station
%
% DEPENDENCIES:
%   - t_tide.m (see Deltares' OpenEarthTools)
%
% This routine was developed as part of the TKI ShorelineS

%% Version
%
% ORIGINAL AUTHOR
% Arvid Dujardin; AnteaGroup Belgium
%
% DATE
% Oct-2023; Feb-2024
%
% $Id: $
% $Date: $
% $Author: $
% $Revision: $
% $HeadURL: $

%% check input

% check dimensions xystat
[N,M]=size(xystat);
if M~=2
    error('xystat should be a [Nx2] matrix')
elseif N<2
    error('xystat should list the coordinates of at least 2 stations')
end

% check dimensions cZWL
if ~iscell(cZWL)
    error('zwl should be a cell-array')
elseif length(cZWL)~=N
    error(['cZWL should have as many cells as waterlevel stations',...
           'defined in xystat'])
else
    [~,nColumns]=cellfun(@size,cZWL);
    if unique(nColumns)~=2
        error('cZWL should only contain [nTx2] matrices')
    end
end

% check dimensions vmean
if numel(vmean)==1
    vmean=repmat(vmean,[N,1]);
elseif length(vmean)~=N
    error(['vmean should be a vector, stating the tidally averaged',...
           '(residual) alongshore velocity for every station'])
end

% check dimensions xsdst
if numel(xsdst)==1
    xsdst=repmat(xsdst,[N,1]);
elseif length(xsdst)~=N
    error(['xsdst should be a vector, stating the distance between the '...
           'station and the shoreline'])
end

% check dimensions dps
if numel(dps)==1
    dps=repmat(dps,[N,1]);
elseif length(dps)~=N
    error('dps should be a vector, stating the depth in every station')
end

%% Handling varargin keyword-value pairs

% set defaults
latSwitch = 0; % false
ssSwitch = '1D_analytical';
Cf = 9.81/65^2;           
hmin = 0.1;
tidefile = 'tidefile.txt';

% if varargin is given, overwrite defaults
while ~isempty(varargin)
    if ischar(varargin{1})
        % keyword-value pairs
        switch varargin{1}
            case 'latitude'
                latSwitch = 1;
                lat=varargin{2};
                % check dimensions lat
                if numel(lat)==1
                    lat=repmat(lat,[N,1]);
                elseif length(lat)~=N
                    error(['lat should be a vector, stating the ',...
                           'latitude for every station'])
                end
            case 'ssSwitch'
                ssSwitch=varargin{2};
                if ~strcmp(ssSwitch,'1D_analytical') && ...
                   ~strcmp(ssSwitch,'harmonic')
                    error(['value for ssSwitch can only be ',...
                           '''1D_analytical'' or ''harmonic''.']);                 
                end
            case 'Cf'
                Cf=varargin{2};
            case 'hmin'
                hmin=varargin{2};
            case 'tidefile'
                tidefile=varargin{2};
            otherwise
                error(['Unknown keyword: ',varargin{1}]);
        end
        varargin([1 2])=[];
    else
        error('varargin should be given as keyword-value pairs')
    end
end

%% Harmonic analysis

disp(['Performing harmonic analysis on ',num2str(N),' waterlevel stations'])
cTide=cell(N,1);
dt=NaN(N,1);
for n=1:N
    time=cZWL{n}(:,1);
    % check if timeseries is equidistant with single-precision tolerance
    % (uniquetol scales the tol input based on the magnitude of the data)
    dt_tmp=uniquetol(diff(time),1e-6); % single-precision tolerance
    if numel(dt_tmp)~=1
        error(['waterlevel timeseries at station ',num2str(n),...
               ' should be equidistant in time'])
    else
        dt(n)=dt_tmp;
    end
    if ~latSwitch % raw constituent phases at the central time
        cTide{n}=t_tide(cZWL{n}(:,2),'interval',dt(n)*24,...
                        'error','wboot','output','none','diary','none');
    else % t_tide with nodal correction
        cTide{n}=t_tide(cZWL{n}(:,2),'interval',dt(n)*24,...
                        'latitude',lat(n),'start time',datevec(time(1)),...
                        'error','wboot','output','none','diary','none');
    end
end

%% Preparing output

disp('Preparing ShorelineS tidefile.txt')
% Dano's time2har.m produces eta, phasewl, periodwl and offsetwl
% comparison with t_tide.m:
% => eta      = tidecon(i,1)
%    phasewl  = tidecon(i,3)
%    periodwl = 1/freq
%    offsetwl = z0

% output in tidefile:
% [xstat ystat etaM2 etaM4 detadsM2 detadsM4 phiM2 phiM4 kM2 kM4 ss]
M2i=NaN(N,1);
M4i=NaN(N,1);
eta=NaN(N,2);
phi=NaN(N,2);
offsetwl=NaN(N,1);
for n=1:N
    M2i(n)=find(strcmp(cellstr(cTide{n}.name),'M2'));
    M4i(n)=find(strcmp(cellstr(cTide{n}.name),'M4'));
    eta(n,:)=cTide{n}.tidecon([M2i(n) M4i(n)],1)';
    phi(n,:)=cTide{n}.tidecon([M2i(n) M4i(n)],3)';
    offsetwl(n)=cTide{n}.z0;
end
% For the first and the last station deltaS, deltaPhi and deltaEta/deltaS 
% are based on the differences between the current and the next/previous
% station. For all other stations deltaS, deltaPhi and deltaEta/deltaS are 
% based on the differences between the previous and the next station 
xstat=xystat(:,1);
ystat=xystat(:,2);
ds=NaN(N,1);
dphi=NaN(N,2);
k=NaN(N,2);
detads=NaN(N,2);
surfslopem=NaN(N,1);
surfslopec=NaN(N,1);
for n=1:N
    im=max(n-1,1);
    ip=min(n+1,length(xstat));
    ds(n)=hypot(xstat(ip)-xstat(im),ystat(ip)-ystat(im));
    dphi(n,:)=(mod(phi(ip,:)-phi(im,:)+180,360)-180)*pi/180;
    k(n,:)=dphi(n,:)/ds(n);
    detads(n,:)=(eta(ip,:)-eta(im,:))/ds(n);
    surfslopem(n)=(offsetwl(ip)-offsetwl(im))/ds(n);
    Ttide=1/cTide{n}.freq(M2i(n)); % duration of the M2 tidal cycle in hours
    nT=ceil(Ttide/(dt(n)*24));  % #output timesteps per M2 tidal cycle
    if strcmp(ssSwitch,'1D_analytical')
        % get analytical solution for longshore surface slope instead of observed one
        [~,~,surfslopec(n)]=tide_1d_ana_anycomp(eta(n,:),detads(n,:),...
                                                phi(n,:),vmean(n),...
                                                surfslopem(n),Ttide*60*60,...
                                                nT,k(n,:),Cf,hmin,...
                                                xsdst(n),dps(n));
    else
        % use longshore surface slope from harmonic analysis
        surfslopec(n)=surfslopem(n);
    end
end

%% output txt file tidal forcing for ShorelineS

disp('Writing ShorelineS tidefile')
out=[xystat,eta,detads,phi,k,surfslopec];
fid=fopen(tidefile,'w');
fprintf(fid,'%8s %8s %7s %7s %11s %11s %5s %5s %11s %11s %11s\n',...
        '%  xstat','ystat','etaM2','etaM4','detadsM2','detadsM4','phiM2','phiM4','kM2','kM4','ss');
for n=1:N
    fprintf(fid,'%8.0f %8.0f %7.3f %7.3f %11.3d %11.3d %5.1f %5.1f %11.3d %11.3d %11.3d',...
            out(n,:));
    fprintf(fid,'\n');
end
fclose(fid);

end % function
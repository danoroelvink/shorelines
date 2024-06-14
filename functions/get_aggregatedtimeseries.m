function [D]=get_aggregatedtimeseries(D,TIME)
% function [D]=get_aggregatedtimeseries(D,TIME)
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
    
    
    %% AGGREGATE TIME-SERIES DATA
    % method of speeding up the simulation by computing average conditions for time step periods
    if ~isempty(D(1).timenum) && TIME.tc==0
        fldvc = {'Dir' 'PHIequi' 'PHIf'};
        fldsc = {'Hs' 'Tp' 'c1' 'c2' 'h0' 'fshape' 'QSoffset'};
        sfac = [2.0 1.0 1.0 1.0 1.0 1.0 1.0]; %scaling factors for scalars  
        timestep = TIME.dt*365;
        
        for kk=1:length(D)
            DT=diff(D(kk).timenum);             
            DT0=median(DT);
            factor=max(ceil(timestep/DT0),1);
            DT1=factor*DT0;
            ttend=min(ceil((D(kk).timenum(end)-D(kk).timenum(1))/DT1)+1,length(D(kk).timenum));
            
            if factor>2 
                Dtmp=struct;
                ttrange0=1;
                for tt=1:ttend
                    tstart=D(kk).timenum(1)+(tt-1)*DT1;
                    tend=D(kk).timenum(1)+tt*DT1;
                    if tt==ttend
                        tend=inf;
                    end
                    
                    if max(DT)==min(DT)
                        % for equidistant time-series
                        ttrange=unique(min(max((tt-1)*factor+[-floor((factor-1)/2):ceil((factor-1)/2)],1),length(D(kk).timenum)));
                    else
                        % for non-equidistant time-series
                        % aggregate scalars over period 'ttrange']
                        ttrange=find(D(kk).timenum>=tstart & D(kk).timenum<tend);
                    end  
                    if isempty(ttrange)
                        ttrange=ttrange0(end);
                    end
                    ttrange0=ttrange;
                    
                    Dtmp.timenum(tt,1) = D(kk).timenum(1)+DT1*(tt-1);

                    % aggregate scalars over period 'ttrange'
                    for jj=1:length(fldsc)
                        if isfield(D(kk),fldsc{jj})
                            Dtmp.(fldsc{jj})(tt,1) = (mean(D(kk).(fldsc{jj})(ttrange).^sfac(jj))).^(1./sfac(jj));
                        end
                    end
                    % aggregate vectors over period 'ttrange'
                    for jj=1:length(fldvc)
                        if isfield(D(kk),fldvc{jj})
                            HSfac = sfac(1);
                            sdir0=sind(D(kk).(fldvc{jj})(ttrange)); sdir0 = sdir0(:);
                            cdir0=cosd(D(kk).(fldvc{jj})(ttrange)); cdir0 = cdir0(:);
                            sdir1=mean(sdir0.*D(kk).Hs(ttrange).^HSfac ./ sum(D(kk).Hs(ttrange).^HSfac) );
                            cdir1=mean(cdir0.*D(kk).Hs(ttrange).^HSfac ./ sum(D(kk).Hs(ttrange).^HSfac) );
                            Dtmp.(fldvc{jj})(tt,1) = mod(atan2d(sdir1,cdir1),360);
                        end
                    end
                end
                
                D(kk).timenum      = Dtmp.timenum;
                for jj=1:length(fldsc)
                    if isfield(D(kk),fldsc{jj})
                        D(kk).(fldsc{jj}) = Dtmp.(fldsc{jj});
                    end
                end
                for jj=1:length(fldvc)
                    if isfield(D(kk),fldvc{jj})
                        D(kk).(fldvc{jj}) = Dtmp.(fldvc{jj});
                    end
                end
            end
        end
    end
end

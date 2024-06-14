function [ Hg, phiwg,iwtw,phiwg_wr,Hg_wr  ] = get_interpolated_wavefield_Hs_dir_Tp( xg,yg,Hg_all,phiwg_all,Hso,phiw0,Tp0,Hstab,dirtab,Tptab)
% function [ Hg, phiwg,iwtw,phiwg_wr,Hg_wr  ] = get_interpolated_wavefield_dir_Tp( xg,yg,Hg_all,phiwg_all,Hso,phiw0,Tp0,dirtab,Tptab)
% 
% GET_INTERPOLATED_WAVEFIELD - Interpolates a wave field (Hs, dir) based on
% series of wave fields for different offshore wave heights and directions
% find out the wave conditions in the transformation matrix;
% convention: nHs times nphiw conditions
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

    iwtw=0;
    %phiw0=phiw0*180/pi;
    
    for i=2:length(dirtab)
        if dirtab(i)<dirtab(i-1)
            dirtab(i)=dirtab(i)+360;
        end
    end
    if phiw0<dirtab(1)
        phiw0=phiw0+360;
        if phiw0>dirtab(end)
            disp('phiw0 outside range of dirtab!')
            phiw0=0;
            Hso=0;
            iwtw=1;
        end
    end
    if phiw0>dirtab(end)
        phiw0=phiw0-360;
        if phiw0<dirtab(1)
            disp('phiw0 outside range of dirtab!')
            phiw0=0;
            Hso=0;
            iwtw=1;
        end
    end
    Hso=max(min(Hso,Hstab(end)),Hstab(1));
    phiw0=max(min(phiw0,dirtab(end)),dirtab(1));
    Tp0=max(min(Tp0,Tptab(end)),Tptab(1));
    nHs=length(Hstab);
    ndir=length(dirtab);
    nTp=length(Tptab);
    numHs=[1:nHs];
    numdir=[1:ndir];
    numTp=[1:nTp];
    indHs=interp1(Hstab,numHs,Hso);
    inddir=interp1(dirtab,numdir,phiw0);
    indTp=interp1(Tptab,numTp,Tp0);
    ind(1,1)=min(floor(indHs),nHs-1);
    ind(1,2)=min(floor(inddir),ndir-1);
    ind(1,3)=min(floor(indTp),nTp-1);
    ind(2,:)=ind(1,:)+1;
    w(1,1)=ind(2,1)-indHs;
    w(1,2)=ind(2,2)-inddir;
    w(1,3)=ind(2,3)-indTp;
    w(2,:)=1-w(1,:);
    
    Hg=zeros(size(xg));
    cosphiwg=zeros(size(xg));
    sinphiwg=zeros(size(xg));
    
    for k=1:2
        for j=1:2
            for i=1:2
                i1=ind(i,1);
                j2=ind(j,2);
                k3=ind(k,3);
                Hg=Hg+squeeze(Hg_all(:,:,  i1,    j2,    k3)) ...
                                       *w(i,1)*w(j,2)*w(k,3);
                cosphiwg=cosphiwg+squeeze(cosd(phiwg_all(:,:,  i1,    j2,    k3))  ...
                                                .*Hg_all(:,:,  i1,    j2,    k3)) ...
                                                           *w(i,1)*w(j,2)*w(k,3);
                sinphiwg=sinphiwg+squeeze(sind(phiwg_all(:,:,  i1,    j2,    k3))  ...
                                                .*Hg_all(:,:,  i1,    j2,    k3)) ...
                                                           *w(i,1)*w(j,2)*w(k,3);
            end
        end
    end
    
    phiwg=mod(atan2d(sinphiwg,cosphiwg),360);
    
    
    % row=[22,22,24,30,38];
    % col=[40,46,52,60,63];
    %  
    % for iwr=1:5
    %     phiwg_wr(iwr)=phiwg(row(iwr),col(iwr));
    %     Hg_wr(iwr)=Hg(row(iwr),col(iwr));
    %     % plot(xg(1,col(iwr)),yg(row(iwr),1)
    % end
    phiwg_wr=[];
    Hg_wr=[];
end

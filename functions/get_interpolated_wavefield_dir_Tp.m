function [ Hg, phiwg,iwtw,phiwg_wr,Hg_wr  ] = get_interpolated_wavefield_dir_Tp( xg,yg,Hg_all,phiwg_all,Hso,phiw0,Tp0,dirtab,Tptab)
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
    phiw0=max(min(phiw0,dirtab(end)),dirtab(1));
    nTp=length(Tptab);
    ndir=length(dirtab);
    numTp=[1:nTp];
    numdir=[1:ndir];
    % indTp=interp1(Tptab,numTp,Tp0,'extrap');
    % inddir=interp1(dirtab,numdir,phiw0,'extrap');
    indTp=interp1(Tptab,numTp,Tp0);
    inddir=interp1(dirtab,numdir,phiw0);
    
    if any(isnan(indTp))||any(isnan(inddir))
        error('get_interpolated_wavefield_dir_Tp :: Requested Tp or direction falls out of range of wave lookup table. Change wavecon table, or rerun wave model.'); 
    end
    
    iT1=min(floor(indTp),nTp-1);
    iT2=iT1+1;
    wT1=iT2-indTp;
    wT2=1-wT1;
    id1=min(floor(inddir),ndir-1);
    id2=id1+1;
    wd1=id2-inddir;
    wd2=1-wd1;
    
    w1=wT1*wd1;
    w2=wT2*wd1;
    w3=wT1*wd2;
    w4=wT2*wd2;
    
    Hg=squeeze(Hg_all(:,:,id1,iT1)*w1+Hg_all(:,:,id1,iT2)*w2 ...
        +Hg_all(:,:,id2,iT1)*w3+Hg_all(:,:,id2,iT2)*w4);
    %phiwg=squeeze(phiwg_all(:,:,id1,iT1)*w1+phiwg_all(:,:,id1,iT2)*w2 ...
    %    +phiwg_all(:,:,id2,iT1)*w3+phiwg_all(:,:,id2,iT2)*w4);
    cosphiwg=squeeze(cosd(phiwg_all(:,:,id1,iT1)).*Hg_all(:,:,id1,iT1)*w1 ...
                    +cosd(phiwg_all(:,:,id1,iT2)).*Hg_all(:,:,id1,iT2)*w2 ...
                    +cosd(phiwg_all(:,:,id2,iT1)).*Hg_all(:,:,id2,iT1)*w3 ...
                    +cosd(phiwg_all(:,:,id2,iT2)).*Hg_all(:,:,id2,iT2)*w4);
    sinphiwg=squeeze(sind(phiwg_all(:,:,id1,iT1)).*Hg_all(:,:,id1,iT1)*w1 ...
                    +sind(phiwg_all(:,:,id1,iT2)).*Hg_all(:,:,id1,iT2)*w2 ...
                    +sind(phiwg_all(:,:,id2,iT1)).*Hg_all(:,:,id2,iT1)*w3 ...
                    +sind(phiwg_all(:,:,id2,iT2)).*Hg_all(:,:,id2,iT2)*w4);
    phiwg=atan2d(sinphiwg,cosphiwg);
    
    
    %% scale wave height with Hs0
    Hg=Hg*max(Hso);
    
    row=[22,22,24,30,38];
    col=[40,46,52,60,63];
     
    for iwr=1:5
        phiwg_wr(iwr)=phiwg(row(iwr),col(iwr));
        Hg_wr(iwr)=Hg(row(iwr),col(iwr));
        % plot(xg(1,col(iwr)),yg(row(iwr),1)
    end
    
end

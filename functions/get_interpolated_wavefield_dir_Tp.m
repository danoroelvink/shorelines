function [ Hg, phiwg ] = get_interpolated_wavefield_dir_Tp( xg,yg,Hg_all,phiwg_all,Hs0,phiw0,Tp0,dirtab,Tptab)
%GET_INTERPOLATED_WAVEFIELD - Interpolates a wave field (Hs, dir) based on
% series of wave fields for different offshore wave heights and directions

%% find out the wave conditions in the transformation matrix;
%% convention: nHs times nphiw conditions

for i=2:length(dirtab)
    if dirtab(i)<dirtab(i-1)
        dirtab(i)=dirtab(i)+360;
    end
end
if phiw0<dirtab(1)
    phiw0=phiw0+360;
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
phiwg=squeeze(phiwg_all(:,:,id1,iT1)*w1+phiwg_all(:,:,id1,iT2)*w2 ...
    +phiwg_all(:,:,id2,iT1)*w3+phiwg_all(:,:,id2,iT2)*w4);
%% scale wave height with Hs0
Hg=Hg*max(Hs0);
end


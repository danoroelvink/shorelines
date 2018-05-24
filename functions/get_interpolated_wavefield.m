function [ Hg, phiwg ] = get_interpolated_wavefield( xg,yg,Hg_all,phiwg_all,Hs0,phiw0)
%GET_INTERPOLATED_WAVEFIELD - Interpolates a wave field (Hs, dir) based on
% series of wave fields for different offshore wave heights and directions

%% find out the wave conditions in the transformation matrix;
%% convention: nHs times nphiw conditions
ntot=size(Hg_all,1);
for ii=1:ntot
    if Hg_all(ii,1,1)>Hg_all(ii+1,1,1)
        nHs=ii;
        break
    end
end
ndir=ntot/nHs;
Hstab=squeeze(Hg_all(1:nHs,1,1));  %squeeze:Remove singleton dimensions
phiwtab=squeeze(phiwg_all(1:nHs:ntot,1,1));

for i=2:length(phiwtab)
    if phiwtab(i)<phiwtab(i-1)
        phiwtab(i)=phiwtab(i)+360;
    end
end
Hs0=max(min(Hs0,Hstab(end)),Hstab(1));
if phiw0<phiwtab(1)
    phiw0=phiw0+360;
end
phiw0=max(min(phiw0,phiwtab(end)),phiwtab(1));
numHs=[1:nHs];
numdir=[1:ndir];
indHs=interp1(Hstab,numHs,Hs0); %%
inddir=interp1(phiwtab,numdir,phiw0);
iH1=min(floor(indHs),nHs-1);
iH2=iH1+1;
wH1=iH2-indHs;
wH2=1-wH1;
id1=min(floor(inddir),ndir-1);
id2=id1+1;
wd1=id2-inddir;
wd2=1-wd1;

i1=iH1+(id1-1)*nHs;
i2=iH2+(id1-1)*nHs;
i3=iH1+(id1  )*nHs;
i4=iH2+(id1  )*nHs;
w1=wH1*wd1;
w2=wH2*wd1;
w3=wH1*wd2;
w4=wH2*wd2;

Hg=squeeze(Hg_all(i1,:,:)*w1+Hg_all(i2,:,:)*w2+Hg_all(i3,:,:)*w3+Hg_all(i4,:,:)*w4);
phiwg=squeeze(phiwg_all(i1,:,:)*w1+phiwg_all(i2,:,:)*w2+phiwg_all(i3,:,:)*w3+phiwg_all(i4,:,:)*w4);
end


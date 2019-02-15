clear all;close all
W0.gridfile='RiaFormosa\Ria_rec_grid_UTM_rot.grd';
W0.depfile='RiaFormosa\Ria_rec_dep_UTM_rot.dep';
W0.wavefile='RiaFormosa\Ria_Formosa';
W0.ntheta=18;           % Number of directional bins
W0.thetamin=90;         % In Nautical convention: 0 i FROM north
W0.thetamax=270;        % In Nautical convention: 0 i FROM north
W0.dirmin=110;
W0.dirstep=20;
W0.dirmax=250;
W0.Tmin=6;
W0.Tstep=2;
W0.Tmax=14;
W0.anim=0;
W=waverefrac_implicit(W0);
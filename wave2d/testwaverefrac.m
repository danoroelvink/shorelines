clear all;close all
W0.gridfile='rotated.grd';
W0.depfile='smooth_rotated.dep';
W.wavefile='Damietta_waves.mat';
W=waverefrac_implicit(W0);
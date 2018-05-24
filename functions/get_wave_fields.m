function [ xg,yg,Hg,dirg,dirtab,Tptab ] = get_wave_fields_from_mat( fname )
% GET_WAVE_FIELDS
% routine to generate a series of wave fields for given wave
% height, wave directions and wave periods. The conditions are given  
% in jonswaptable.txt. Results computed by XBeach are stored in 
% xboutput.nc and read into the matrices xg,yg for the grid and
% Hg_all and phiwg_all; these represent a transformation matrix
% from which the wave field for an arbitrary wave condition can be
% interpolated.

%% For now assume that the wave model has run and produced nHs by nphi 
%% wave conditions
load(fname)

end


function [COAST] = get_active_profile(COAST)
%UNTITLED2 Summary of this function goes here
%INPUT:
%   COAST          Structure with data of the ShorelineS model, of which is used
%    .x          x-coordinate of coastal segment [m]
%    .y          y-coordinate of coastal segment [m]
%    .n          number of grid cells of coastal segment
%    .h0       Active profile height input file (scalar or table with 3 columns, h0x, h0y and h0)
%
% OUTPUT:
%   COAST
%      .h0       Active profile height [m]
%      .h0_x     Active profile x postition
%      .h0_y     Active profile y postition
%

    if ischar(COAST.h0input)   % input from file
        % read table with active profile heights at point locations (with 3 columns, h0x, h0y and h0)
        % INPUT EXAMPLE :   S.d='ActiveProfile.txt', 3 column text
        % file
        h0input=load(COAST.h0input);
        h0_x=COAST.h0input(:,1)';
        h0_y=COAST.h0input(:,2)';
        h0=COAST.h0input(:,3)';
        
        % find the right alongshore location for each of the active profile
        % interpolate h0         
        var1=struct;
        var2=struct;
        var2.h0_mc=h0;
        [~,var2i,~]=get_interpolation_on_grid('weighted_distance',COAST.x,COAST.y,h0_x,h0_y,var1,var2);        
        h0=var2i.h0_mc;  
        h0_x=COAST.x;
        h0_y=COAST.y;
        
    elseif isscalar(COAST.h0input)
        % Active profile at fixed height for the whole grid
        % INPUT EXAMPLE :   COAST.PHIf0=10;
        h0_x=COAST.x;
        h0_y=COAST.y;
        h0=repmat(COAST.h0input,[1,COAST.n]);

    elseif isnumeric(COAST.h0input)  % array from input structure
        h0_x=COAST.h0input(:,1)';
        h0_y=COAST.h0input(:,2)';
        h0=COAST.h0input(:,3)';

        % find the right alongshore location for each of the active profile
        % interpolate h0         
        var1=struct;
        var2=struct;
        var2.h0_mc=h0;
        [~,var2i,~]=get_interpolation_on_grid('weighted_distance',COAST.x,COAST.y,h0_x,h0_y,var1,var2);        
        h0=var2i.h0_mc;  
        h0_x=COAST.x;
        h0_y=COAST.y;
        
    else     
        error('get_active_profile::Invalid active profile input');    
    end
    
    % Output to structure
    COAST.h0=h0;
    COAST.h0_x=h0_x;
    COAST.h0_y=h0_y;
   
end
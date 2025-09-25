function [CHANNEL]=prepare_channel(S)
% function [CHANNEL]=prepare_channel(S)
%
% Initialize the data-structure for the channel migration function ('CHANNEL').
%
% INPUT: 
%    S
%         .channel          : switch for migrating inlet (0/1)
%         .channelwidth     : target channel width [m]
%         .channelfac       : adaptation factor, scales the response of the channel [-]
%         .channeldischrate : discharge rate
%         .channeldischr    : discharge rate of river
%         .ldbchannel       : file with initial channel axis x and y coordinates, or directly a [Nx2] matrix (option 1)
%         .xrmc             : x-coordinates of rivers [m] (option 2)
%         .yrmc             : y-coordinates of rivers [m] (option 2)
%         .xyoffset         : offset of x and y axis used for plotting [1x2]
%         
% OUTPUT:
%    CHANNEL
%         .xrmc             : x-coordinates of rivers [m]
%         .yrmc             : y-coordinates of rivers [m]
%         .width            : width of the channel [m]
%         .fac              : factor for the reshaping of the channel, scales the response of the channel [-]
%         .disch_rate       : discahrge rate
%         .disch_R          : discharge rate of river
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
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses>
%   --------------------------------------------------------------------

    fprintf('  Prepare channel \n');
    
    %% Migrating channel
    CHANNEL=struct;
    CHANNEL.used=S.channel;
    CHANNEL.xrmc=[];
    CHANNEL.yrmc=[];
    CHANNEL.width=[];
    CHANNEL.fac=[];
    CHANNEL.disch_rate=[];
    CHANNEL.disch_R=[];
    
    if S.channel
        % SET THE X,Y LOCATION
        if strcmpi(S.ldbchannel,'manual') || strcmpi(S.ldbchannel,'interactive') 
            figure(11);
            xl=xlim;yl=ylim;
            htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add channel axis (LMB); Next channel (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
            [xrmc,yrmc]=select_multi_polygon('k');
            CHANNEL.xrmc=xrmc;
            CHANNEL.yrmc=yrmc;
            set(htxt2,'Visible','off');
            save('rivers.mat','xrmc','yrmc');
        elseif ~isempty(S.xrmc)
            % direct input via a keyword
            CHANNEL.xrmc=S.xrmc;
            CHANNEL.yrmc=S.yrmc;
        elseif ~isempty(S.ldbchannel)
            xy_channel=load(S.ldbchannel);
            if isstruct(xy_channel)
                % not-preferred option using a matlab structure
                CHANNEL.xrmc=xy_channel.xrmc(:)'-S.xyoffset(1);
                CHANNEL.yrmc=xy_channel.yrmc(:)'-S.xyoffset(2);
            else
                % preferred option using [Nx2] matrix for xr and yr
                if size(xy_channel,2)>2
                    xy_channel=xy_channel';
                end
                CHANNEL.xrmc=xy_channel(:,1)'-S.xyoffset(1);
                CHANNEL.yrmc=xy_channel(:,2)'-S.xyoffset(2);
            end
        end
        
        % SET THE WIDTH
        % using predefined values for channel width
        nriv=sum(isnan(CHANNEL.xrmc))+1;
        CHANNEL.width=S.channelwidth;
        CHANNEL.fac=S.channelfac;
        CHANNEL.disch_rate=S.channeldischrate;
        CHANNEL.disch_R=S.channeldischr;
        
        % make sure length is equal to 'nriv'
        if length(CHANNEL.width)<nriv
        CHANNEL.width=repmat(S.channelwidth,[nriv,1]);
        end
        if length(CHANNEL.fac)<nriv
        CHANNEL.fac=repmat(S.channelfac,[nriv,1]);
        end
        if length(CHANNEL.disch_rate)<nriv
        CHANNEL.disch_rate=repmat(S.channeldischrate,[nriv,1]);
        end
        if length(CHANNEL.disch_R)<nriv
        CHANNEL.disch_R=repmat(S.channeldischr,[nriv,1]);
        end
        
        % not-preferred option using a matlab structure (to be removed)
        if ~isempty(S.ldbchannel)
            matname=S.ldbchannel;
            if strcmpi(matname(end-3:end),'.mat')
                channels=struct;
                load(matname);
                n_chan=length(channels);
                for ichan=1:n_chan
                    CHANNEL.width(ichan)= channels(ichan).channelwidth;
                    CHANNEL.fac(ichan)  = channels(ichan).channelfac;
                end
            end
        end
        
        %% write logfile
        % struct2log(CHANNEL,'CHANNEL','a');
    end
end

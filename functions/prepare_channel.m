function [CHANNEL]=prepare_channel(S)
% [xr_mc,yr_mc,channel_width,channel_fac,x_flood,y_flood,flood_deficit,fcell_area]=prepare_channel(S)
%
% INPUT:
%    S
%         .channel        :  
%         .xr_mc          :  
%         .yr_mc          :  
%         .LDBchannel     :  
%         .XYoffset       :  
%         
% OUTPUT:
%    CHANNEL
%         .xr_mc
%         .yr_mc
%         .width
%         .fac
%         .disch_rate
%         .disch_R
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

    fprintf('  Prepare channel \n');
    
    %% Migrating channel
    CHANNEL=struct;
    CHANNEL.used=S.channel;
    CHANNEL.xr_mc=[];
    CHANNEL.yr_mc=[];
    CHANNEL.width=[];
    CHANNEL.fac=[];
    CHANNEL.disch_rate=[];
    CHANNEL.disch_R=[];

    if S.channel
        % SET THE X,Y LOCATION
        if ~isempty(S.xr_mc)
            CHANNEL.xr_mc=S.xr_mc;
            CHANNEL.yr_mc=S.yr_mc;
        elseif ~isempty(S.LDBchannel)
            xy_channel=load(S.LDBchannel);
            try
               CHANNEL.xr_mc=xy_channel(:,1)'-S.XYoffset(1);
               CHANNEL.yr_mc=xy_channel(:,2)'-S.XYoffset(2);
            catch
                CHANNEL.xr_mc=xy_channel.xr_mc;
                CHANNEL.yr_mc=xy_channel.yr_mc;
            end
        else
            figure(11);
            plot(S.x_mc,S.y_mc,'linewidth',2)
            hold on
            xl=xlim;yl=ylim;
            htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add channel axis (LMB); Next channel (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
            [xr_mc,yr_mc]=select_multi_polygon('k');
            CHANNEL.xr_mc=xr_mc;
            CHANNEL.yr_mc=yr_mc;
            set(htxt2,'Visible','off');
            save('rivers.mat','xr_mc','yr_mc');
        end
        
        % SET THE WIDTH
        if ~isempty(S.LDBchannel)
            matname=S.LDBchannel;
            matname(end-2:end)='mat';
            if exist(matname)==2
                channels=struct;
                load(matname);
                n_chan=length(channels);
                for ichan=1:n_chan
                    CHANNEL.width(ichan)= channels(ichan).channel_width;
                    CHANNEL.fac(ichan)  = channels(ichan).channel_fac;
                end
            else
                CHANNEL.width=S.channel_width;
                CHANNEL.fac=S.channel_fac;
                CHANNEL.disch_rate=S.channel_disch_rate;
                CHANNEL.disch_R=S.channel_disch_R;
            end
        else
            nriv=sum(isnan(CHANNEL.xr_mc))+1;
            CHANNEL.width=repmat(S.channel_width,[nriv,1]);
            CHANNEL.fac=repmat(S.channel_fac,[nriv,1]);
            CHANNEL.disch_rate=repmat(S.channel_disch_rate,[nriv,1]);
            CHANNEL.disch_R=repmat(S.channel_disch_R,[nriv,1]);
        end
        %% write logfile
        % struct2log(CHANNEL,'CHANNEL','a');

    end
end

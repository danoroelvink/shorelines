function [COAST,TRANSP,WAVE]=get_spikesremoved(COAST,TRANSP,WAVE)
% function [COAST,TRANSP,WAVE]=get_spikesremoved(COAST,TRANSP,WAVE);
%
% removes spiky grid cells in the x,y coastline variable
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
    if ~isempty(COAST.Wberm)
        return;
    end
    
    ds0=COAST.ds0;
    if ~isscalar(ds0)
        ds0=min(ds0(:,3));
    end
    thr_fraction_of_gridsize=0.2;

    n_mc=length(find(isnan(COAST.x_mc)))+1;
    for i_mc=n_mc:-1:1
        [ x,y,~,i1,i2 ] = get_one_polygon( COAST.x_mc,COAST.y_mc,i_mc );
        ar=0;
        s=zeros(size(x));
        
        if length(x)>1
            for i=1:length(x)
                if i>=2
                    s(i)=s(i-1)+hypot(x(i)-x(i-1),y(i)-y(i-1));
                end
                if i==1
                    s(1)=0;
                    ds=hypot(x(2)-x(end),y(2)-y(end));
                elseif i>1 && i<length(x)
                    ds(i)=hypot(x(i+1)-x(i-1),y(i+1)-y(i-1));
                elseif i==length(x)               
                    ds(i)=hypot(x(1)-x(end-1),y(1)-y(end-1));
                end
            end

            IDspiky=(ds<thr_fraction_of_gridsize*ds0);
%             x(IDspiky)=nan;
%             y(IDspiky)=nan;
%             COAST.x_mc(i1:i2)=x;
%             COAST.y_mc(i1:i2)=y;
            xnew=x(~IDspiky);
            ynew=y(~IDspiky);
            [COAST.x_mc,COAST.y_mc,i1,i2]=insert_section(xnew,ynew,COAST.x_mc,COAST.y_mc,i_mc);
            
            if isempty(xnew)
                try
                    [COAST.dSds_mc]=insert_section(xnew,COAST.dSds_mc,i_mc);
                    [COAST.PHIc_mc]=insert_section(xnew,COAST.PHIc_mc,i_mc);
                    [COAST.PHIf_mc]=insert_section(xnew,COAST.PHIf_mc,i_mc);
                    [TRANSP.QS_mc]=insert_section(xnew,TRANSP.QS_mc,i_mc);
                    [WAVE.HSo_mc]=insert_section(xnew,WAVE.HSo_mc,i_mc);
                    [WAVE.PHIo_mc]=insert_section(xnew,WAVE.PHIo_mc,i_mc);
                    [WAVE.TP_mc]=insert_section(xnew,WAVE.TP_mc,i_mc);
                    [WAVE.HStdp_mc]=insert_section(xnew,WAVE.HStdp_mc,i_mc);
                    [WAVE.PHItdp_mc]=insert_section(xnew,WAVE.PHItdp_mc,i_mc);
                    [WAVE.dPHItdp_mc]=insert_section(xnew,WAVE.dPHItdp_mc,i_mc);
                    [WAVE.HSbr_mc]=insert_section(xnew,WAVE.HSbr_mc,i_mc);
                    [WAVE.PHIbr_mc]=insert_section(xnew,WAVE.PHIbr_mc,i_mc);
                    [WAVE.dPHIbr_mc]=insert_section(xnew,WAVE.dPHIbr_mc,i_mc);
                    [WAVE.hbr_mc]=insert_section(xnew,WAVE.hbr_mc,i_mc);
                catch
                    fprintf('Warning : Empty xy-section could not be achieved for transport and waves for i_mc=%1.0f.\n',i_mc);
                end
            end
            COAST.n_mc=COAST.n_mc-1;
            
            %% make transport points xq_mc and yq_mc
            [COAST]=get_transportpoints(COAST);
        end
    end
    
    [COAST,di]=cleanup_nans(COAST);
    
end

function plot_debug(use_debug,COAST,WAVE,TRANSP,im3,ip3)
% function plot_debug(use_debug,COAST.x,COAST.y,WAVE.dPHItdp,TRANSP.debug.QS2,TRANSP.debug.QS0,TRANSP.debug.QS1,im3,ip3)
%
% INPUT :
%    use_debug  : switch for making the debug plot
%    COAST.x            : x-coordinates of coastline section [Nx1]
%    COAST.y            : y-coordinates of coastline section [Nx1]
%    WAVE.dPHItdp       : relative angle of incidence of the waves [Nx1]
%    TRANSP.debug.QS2   : transport rates : transport computation, upwind correction and shadows [Nx1]
%    TRANSP.debug.QS0   : transport rates : transport computation [Nx1]
%    TRANSP.debug.QS1   : transport rates : transport computation and upwind correction [Nx1]
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

    if use_debug==1
        %if length(x)>100
        fc9=gcf;
        figure(112);
        clf;plot_figureproperties(gcf,1220,920,30);
        
        x2=(COAST.x(1:end-1)+COAST.x(2:end))/2;
        y2=(COAST.y(1:end-1)+COAST.y(2:end))/2;
        
        %% dPHI
        subplot(2,2,1);
        plot(COAST.x,COAST.y);hold on;
        title('WAVE.dPHItdp [?]');
        scatter(x2,y2,10,WAVE.dPHItdp);
        set(gca,'clim',[0,360]);
        colorbar;
        try;plot(COAST.x(im3(2:end,1)),COAST.y(im3(2:end,1)),'ko');end
        %xlim(1.0e+04 * [ 7.1911    7.3864 ]);
        %ylim(1.0e+05 * [ 4.5184    4.5359 ]);
        
        %% TRANSP.debug.QS2
        subplot(2,2,2);
        plot(COAST.x,COAST.y);hold on;
        title('TRANSP.debug.QS2');
        scatter(x2,y2,10,TRANSP.debug.QS2);
        colorbar;
        try;plot(COAST.x(im3(2:end,1)),COAST.y(im3(2:end,1)),'ks');end;
        %xlim(1.0e+04 * [ 7.1911    7.3864 ]);
        %ylim(1.0e+05 * [ 4.5184    4.5359 ]);
        %TRANSP.debug.QS2(IDshadow1)=0.;
        %TRANSP.debug.QS2(IDshadow2)=0.;
        
        %% TRANSP.debug.QS0
        subplot(2,2,3);
        plot(COAST.x,COAST.y);hold on;
        title('TRANSP.debug.QS0');
        scatter(x2,y2,10,TRANSP.debug.QS0);
        colorbar;
        try;plot(x2(im3),y2(im3),'ks');end
        %xlim(1.0e+04 * [ 7.1911    7.3864 ]);
        %ylim(1.0e+05 * [ 4.5184    4.5359 ]);
        
        %% TRANSP.debug.QS1
        subplot(2,2,4);
        plot(COAST.x,COAST.y);hold on;
        title('TRANSP.debug.QS1');
        scatter(x2,y2,10,TRANSP.debug.QS1);
        colorbar;
        try;plot(x2(ip3),y2(ip3),'k+');end
        %set(gca,'clim',[0,1]);colorbar;
        %xlim(1.0e+04 * [ 7.1911    7.3864 ]);
        %ylim(1.0e+05 * [ 4.5184    4.5359 ]);
        
        figure(fc9);
        %end
    elseif use_debug==3
        figure;subplot(4,1,1);plot(TRANSP.QS);grid on;subplot(4,1,2);plot(COAST.yq);grid on;subplot(4,1,3);plot(WAVE.HStdp);hold on;plot(WAVE.HSbr,'r');grid on;subplot(4,1,4);plot(WAVE.PHItdp);hold on;plot(WAVE.PHIbr,'r');grid on
    elseif use_debug==4
        hold on;
        hsc=scatter(COAST.xq_mc,COAST.yq_mc,5,TRANSP.QS_mc);
        for kk=1:length(TRANSP.QS_mc)
            hpp=text(COAST.xq_mc(kk),COAST.yq_mc(kk),num2str(TRANSP.QS_mc(kk),'%1.0f'));
        end
    elseif use_debug==5
        % test symmetry
        %figure(88);clf;plot(TRANSP.QS_mc);hold on;plot(fliplr([1:length(TRANSP.QS_mc)]),-TRANSP.QS_mc,'r--')
        figure(88);clf;subplot(3,1,1);plot(TRANSP.QS_mc);hold on;plot(fliplr([1:length(TRANSP.QS_mc)]),-TRANSP.QS_mc,'r--');xi=[1:length(TRANSP.QS_mc)];plot(xi(TRANSP.shadowS_mc),TRANSP.QS_mc(TRANSP.shadowS_mc),'ko');subplot(3,1,2);plot(270-WAVE.PHIbr_mc);hold on;plot(fliplr([1:length(WAVE.PHIbr_mc)]),-(270-WAVE.PHIbr_mc),'r--');subplot(3,1,3);plot(COAST.x_mc,COAST.y_mc);hold on;plot(COAST.x_mc,-COAST.y_mc,'r--');plot(COAST.xq_mc,COAST.yq_mc,'g*');plot(COAST.xq_mc(TRANSP.shadowS_mc),COAST.yq_mc(TRANSP.shadowS_mc),'ko')
        %figure(88);clf;subplot(3,1,1);plot(TRANSP.QS);hold on;plot(fliplr([1:length(TRANSP.QS)]),-TRANSP.QS,'r--');xi=[1:length(TRANSP.QS)];plot(xi(TRANSP.shadowS),TRANSP.QS(TRANSP.shadowS),'ko');subplot(3,1,2);plot(mod(270-WAVE.PHIbr,360));hold on;plot(fliplr([1:length(WAVE.PHIbr)]),mod(-(270-WAVE.PHIbr),360),'r--');subplot(3,1,3);plot(COAST.x_mc,COAST.y_mc);hold on;plot(COAST.x_mc,-COAST.y_mc,'r--');plot(COAST.xq,COAST.yq,'g*');plot(COAST.xq(TRANSP.shadowS),COAST.yq(TRANSP.shadowS),'ko')
        
    end
end
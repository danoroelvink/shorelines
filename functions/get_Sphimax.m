function [WAVE,TRANSP]=get_Sphimax(WAVE,TIDE,TRANSP,STRUC)
% function [WAVE,TRANSP]=get_Sphimax(WAVE,TIDE,TRANSP,STRUC)
% 
% The maximum transport (QSmax) is computed under high-angle wave incidence, 
% as well as the critical wave angle with respect to the coastline at
% the depth-of-closure (.dPHIcrit) and at the point of breaking (.dPHIcritbr). 
% 
% INPUT: 
%    WAVE
%       .HStdp      : wave height at nearshore location (depth-of-closure)
%       .htdp       : water depth at nearshore location (depth-of-closure)
%       .TP         : wave period at nearshore location (depth-of-closure) 
%    TRANSP         : input structure with TRANSPORT data  
%       .gamma      : wave breaking coefficient
% 
% OUTPUT:
%    TRANSP
%       .QSmax      : maximum transport at high-angle wave incidence [m3/yr]
%    WAVE
%       .dPHIcrit   : critical angle at depth-of-closure [°]
%       .dPHIcritbr : critical angle at point of breaking [°]
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


    %% find maximum by fitting a parabola
    eps=1d-4;
    WAVE.dPHIcrit=[];     % at depth-of-closure
    TRANSP.QSmax=[];
    nw=size(WAVE.HStdp,1); % nr of wave conditions
    nq=size(WAVE.HStdp,2); % nr of grid cells of the coast
    QSmax=zeros([nw,nq]);
    waveangle1=35;
    waveangle2=45;
    WAVE.dPHIcrit=nan([nw,nq]);
    WAVE.dPHIcritbr=nan([nw,nq]);
    
    % computation 1
    WAVE1=WAVE;
    WAVE1.dPHItdp=repmat(waveangle1,[nw,nq]);
    WAVE1.HStdp=max(WAVE1.HStdp,eps);
    [WAVE1]=wave_breakingheight(WAVE1,TRANSP);
    WAVE1.HStdp=min(WAVE1.HStdp,WAVE1.ddeep*0.5);
    [TRANSP]=transport(TRANSP,WAVE1,TIDE,STRUC,'QSmax');
    dPHI1=WAVE1.dPHItdp;
    QS1=TRANSP.QSmax; 
    
    % computation 2
    WAVE2=WAVE;
    WAVE2.dPHItdp=repmat(waveangle2,[nw,nq]);
    WAVE2.HStdp=max(WAVE2.HStdp,eps);
    [WAVE2]=wave_breakingheight(WAVE2,TRANSP);
    WAVE2.HStdp=min(WAVE2.HStdp,WAVE2.ddeep*0.5);
    [TRANSP]=transport(TRANSP,WAVE2,TIDE,STRUC,'QSmax');
    dPHI2=WAVE2.dPHItdp;
    QS2=TRANSP.QSmax;
    
    % compute avarage dPHI
    WAVEm=WAVE1;
    dPHIm=.5*(dPHI1+dPHI2);       % at depth-of-closure
    
    % perform iteration for each alongshore grid cell -> map WAVE to WAVE1
    err=repmat(1d3,[nw,nq]);
    iter=0; 
    while ( min(err(:))>eps || iter<=2 ) && iter<10 %err0>err
        iter=iter+1;
        
        % computation 3, 4 etc 
        WAVEm.dPHItdp=dPHIm;
        [WAVEm]=wave_breakingheight(WAVEm,TRANSP);
        [TRANSP]=transport(TRANSP,WAVEm,TIDE,STRUC,'QSmax');
        QSm=TRANSP.QSmax;
        
        % perform iteration for each alongshore grid cell 
        dPHImold=dPHIm;
        in=(err>eps);
        [ivals,jvals]=find(in);
        for i=1:length(ivals)
            i1=ivals(i);
            j1=jvals(i);
            % compute better estimate for dPHIm using computed dPHI's and QS's
            A=[dPHI1(i1,j1)^2,dPHI1(i1,j1),1;dPHI2(i1,j1)^2,dPHI2(i1,j1),1;dPHIm(i1,j1)^2,dPHIm(i1,j1),1];              % at depth-of-closure
            B=[QS1(i1,j1);QS2(i1,j1);QSm(i1,j1)];
            warning off
            a=A\B;
            warning on
            dPHIm(i1,j1)=-a(2)/(2*a(1));     % at depth-of-closure
        end
        
        % renew the dPHI1/dPHI2 and QS1/QS2 for the next iteration cycle
        dPHI2(QS1>QS2)=dPHImold(QS1>QS2);
        dPHI1(QS1<=QS2)=dPHImold(QS1<=QS2);
        QS2(QS1>QS2)=QSm(QS1>QS2);
        QS1(QS1<=QS2)=QSm(QS1<=QS2);
        
        % limit dPHIm
        dPHIm(dPHIm>179.99)=179.99;
        
        % get error value
        err=abs(dPHIm-dPHImold);
        dPHIm(dPHIm==999)=nan;
        %dPHIm(err>10)=nan;
        %QSm(err>10)=nan;
        
        % update the critical angles and maximum transports
        WAVE.dPHIcrit=dPHIm;     % at depth-of-closure
        WAVE.dPHIcritbr=WAVEm.dPHIbr; % at point of breaking
        QSmax=QSm;
    end
    
    % interpolate the critical angles and maximum transports
    for i=1:nw
    WAVE.dPHIcrit(i,:)=interpNANs(WAVE.dPHIcrit(i,:));
    WAVE.dPHIcritbr(i,:)=interpNANs(WAVE.dPHIcritbr(i,:));   
    TRANSP.QSmax(i,:)=interpNANs(QSmax(i,:));
    end

    % special case where no transport could be computed for whole domain (e.g. when wave height is zero)
    WAVE.dPHIcrit(isnan(WAVE.dPHIcrit))=45;
    WAVE.dPHIcritbr(isnan(WAVE.dPHIcritbr))=45;
    TRANSP.QSmax(isnan(TRANSP.QSmax))=0;
    
    % % re-combining the effect of multiple wave conditions to a single value (in case of climate with simultaneous wave conditions)
    % if size(TRANSP.QSmax,1)>1
    %     wghtNR=repmat(WAVE.Prob,[1,size(TRANSP.QSmax,2)]);
    %     HStdp0=max(WAVE.HStdp,1e-6);
    %     wghtHS=HStdp0.^2;
    %     TRANSP.QSmax=sum(TRANSP.QSmax.*wghtNR,1);
    %     sdr = sum(sind(WAVE.dPHIcrit).*wghtNR.*wghtHS,1);
    %     cdr = sum(cosd(WAVE.dPHIcrit).*wghtNR.*wghtHS,1);
    %     WAVE.dPHIcrit=mod(atan2d(sdr,cdr),360);
    %     sdr = sum(sind(WAVE.dPHIcritbr).*wghtNR.*wghtHS,1);
    %     cdr = sum(cosd(WAVE.dPHIcritbr).*wghtNR.*wghtHS,1);
    %     WAVE.dPHIcritbr=mod(atan2d(sdr,cdr),360);
    % end
end

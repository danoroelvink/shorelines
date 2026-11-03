function [WAVE]=wave_breakingheight(WAVE,TRANSP)
% function [WAVE]=wave_breakingheight(WAVE,TRANSP)
% 
% This routine computes the refraction and shoaling from the nearshore point (TDP) to the point of breaking (BR).
% 
% INPUT: 
%    WAVE
%         .dnearshore        : depth at toe of dynamic profile [m]
%         .HStdp             : wave height at nearshore location (or diffracted wave height) [Nx1]
%         .TP                : wave period [s]
%         .dPHItdp           : relative angle of waves with respect to the coastline at nearshore location [°] 
%         .dPHIcrit          : Critical orientation of the coastline with respect to waves at the depth-of-closure (Nx2) [°]
%         .gamma             : depth-induced breaking coefficient [-]
%    TRANSP
%         .suppresshighangle : option to suppress high-angle instabilities by maximizing transport at angles beyond dphicrit (0/1)
%         .trform            : transport formulation (either 'CERC', 'KAMP', 'MILH', 'CERC3', 'VR14')
% 
% OUTPUT:
%    WAVE
%         .HSbr              : breaking wave height (or diffracted wave height) [Nx1]
%         .dPHIbr            : relative angle of waves with respect to the coastline at the point of breaking [°]
%         .dPHItdp           : relative angle of waves with respect to the coastline at nearshore location [°] (updated)
%         .hbr               : depth at point of breaking [m]
%         .cbr               : wave celerity at point of breaking [m/s]
%         .nbr               : deep/shallow water wave number [-]
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

    eps=1d-5;
    hmin=0.1;
    hbrmin=0.05; % minimum depthg of breaking
    WAVE.HSbr=zeros(size(WAVE.dPHItdp));
    WAVE.dPHIbr=zeros(size(WAVE.dPHItdp));
    WAVE.hbr=zeros(size(WAVE.dPHItdp));
    WAVE.cbr=zeros(size(WAVE.dPHItdp));
    WAVE.nbr=zeros(size(WAVE.dPHItdp));
    
    %% Limit high-angle instabilities by forcing maximum transport when the angle is larger than dPHIcrit
    if TRANSP.suppresshighangle==1 
        if isfield(WAVE,'dPHIcrit_mc0')
            dPHIcrit = get_one_polygon( WAVE.dPHIcrit_mc0,WAVE.i_mc);
            if length(dPHIcrit(~isnan(dPHIcrit)))==length(WAVE.dPHItdp)
                WAVE.dPHItdp=sign(WAVE.dPHItdp).*min([abs(WAVE.dPHItdp);dPHIcrit],[],1);        
            elseif length(dPHIcrit(~isnan(dPHIcrit)))>=1
                WAVE.dPHItdp=sign(WAVE.dPHItdp).*min(abs(WAVE.dPHItdp),median(dPHIcrit));
            else
                WAVE.dPHItdp=sign(WAVE.dPHItdp).*min(abs(WAVE.dPHItdp),45);    
            end
        else
            WAVE.dPHItdp=sign(WAVE.dPHItdp).*min(abs(WAVE.dPHItdp),45);        
        end
    end
    
    if ~strcmpi(TRANSP.trform,'RAY') && ~strcmpi(TRANSP.trform,'CERC') && ~strcmpi(TRANSP.trform,'CERC2')
        % Vectorised over all points (and wave conditions): the same secant iteration
        % as the point-by-point loop, applied only to the points that have not
        % converged yet. Identical results, about 30x faster.
        TP=WAVE.TP(:)';
        HStdp=WAVE.HStdp(:)';
        gamma=WAVE.gamma;
        [~,ctdp,~,ntdp]=get_disper(WAVE.dnearshore,TP);
        cosPHI=max(cosd(WAVE.dPHItdp(:)'),eps);
        sinPHI=sind(WAVE.dPHItdp(:)');
        
        % First estimate: Hbr=WAVE.HStdp
        hbr1=max(HStdp./gamma,eps);
        [hbr2,~,err1]=wave_shoalref_vec(hbr1,TP,gamma,HStdp,ctdp,ntdp,sinPHI,cosPHI);
        hbr2b=max(hbr2,eps);
        
        % Second estimate: fill in hbr2
        [hbr3,dPHIbr0,err2,cbr0,nbr0]=wave_shoalref_vec(hbr2b,TP,gamma,HStdp,ctdp,ntdp,sinPHI,cosPHI);
        hbrnew=hbr3;
        
        % Following estimates: inter/extrapolate from last two WAVE.hbr/err pairs
        act=find(~(abs(err2)<eps));
        for iter=3:10
            if isempty(act)
                break
            end
            hn=max(hbr1(act)-err1(act).*(hbr2(act)-hbr1(act))./(err2(act)-err1(act)),eps);
            [~,dPHIa,erra,cbra,nbra]=wave_shoalref_vec(hn,TP(act),gamma,HStdp(act),ctdp(act),ntdp(act),sinPHI(act),cosPHI(act));
            hbrnew(act)=hn;
            dPHIbr0(act)=dPHIa;
            cbr0(act)=cbra;
            nbr0(act)=nbra;
            cont=abs(erra)>eps;
            ac=act(cont);
            hbr1(ac)=hbr2(ac); err1(ac)=err2(ac);
            hbr2(ac)=hn(cont); err2(ac)=erra(cont);
            act=ac;
        end
        
        hbrnew(isnan(hbrnew))=0;
        dPHIbr0(isnan(dPHIbr0))=0;
        sz=size(WAVE.dPHItdp);
        WAVE.HSbr=reshape(hbrnew*gamma,sz);
        WAVE.dPHIbr=reshape(dPHIbr0,sz);
        WAVE.hbr=reshape(hbrnew,sz);
        WAVE.cbr=reshape(cbr0,sz);
        WAVE.nbr=reshape(nbr0,sz);
           
    else 
        % in case TRANSP.trform is 'RAY', 'CERC' or 'CERC2'
        WAVE.HSbr=WAVE.HStdp;
        WAVE.dPHIbr=WAVE.dPHItdp;
        WAVE.hbr=max(WAVE.HStdp./WAVE.gamma,0.1*hbrmin);
        [~,cbr,~,nbr]=get_disper(WAVE.hbr,WAVE.TP);
        WAVE.cbr=cbr;
        WAVE.nbr=nbr;
    end
    
end

function [hbrnew,dPHIbr,err,cbr,nbr]=wave_shoalref_vec(hbr,tper,gamma,hstdp,ctdp,ntdp,sinPHIw,cosPHIw)
% Element-wise (vectorised) version of wave_shoalref with identical branching.
    [~,cbr,~,nbr]=get_disper(hbr,tper);
    dPHIbr=acosd(cosPHIw);   % default: input wave angle (very oblique incidence)
    hstdpbr=hstdp;
    in1=abs(cbr./ctdp.*sinPHIw)<1;
    dPHIbr(in1)=asind(cbr(in1)./ctdp(in1).*sinPHIw(in1));
    cosbr=cosd(dPHIbr);
    in2=in1 & (nbr.*cbr.*cosbr)>0 & abs(cosbr)>1e-3;
    hstdpbr(in2)=hstdp(in2).*sqrt(ntdp(in2).*ctdp(in2).*cosPHIw(in2)./(nbr(in2).*cbr(in2).*cosbr(in2)));
    dPHIbr(in1 & ~in2)=acosd(cosPHIw(in1 & ~in2));
    hbrnew=hstdpbr/gamma;
    err=hbrnew-hbr;
end

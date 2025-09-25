function [COAST,WAVE,TRANSP]=collect_variables(COAST,WAVE,TRANSP,DUNE,MUD,GROYNE,TIME,i_mc)
% function [COAST,WAVE,TRANSP]=collect_variables(COAST,WAVE,TRANSP,DUNE,MUD,GROYNE,TIME,i_mc)
%
% This routine collects the variables of the evaluation of section 'i_mc'
% into the overarching variable that stores data for all sections, with '_mc' suffix. 
% The x_mc and y_mc are temporarily stored in x1_mc and y1_mc, and written 
% back to x_mc and y_mc when the last element has been evaluated.  
% 
% INPUT:
%     COAST   : fields stored 's','ds','x1','y1','xq','yq','xq1','yq1','PHIc','dPHIc','PHIcxy','PHIcxy0','PHIf','h0','cyclic','clockwise'
%     WAVE    : fields stored 'PHIo','PHItdp','PHIbr','HSo','HStdp','HSbr','dPHIo','dPHItdp','dPHIbr','TP','hbr','dPHIcrit','diff'
%     TRANSP  : fields stored 'QS','QSmax','shadowS','shadowS_h','shadowS_hD','shadow','shadow_h','ivals','idrev'
%     DUNE    : fields stored 'qs','qss','ql','qw','R','SWL','wberm','dfelev','dcelev','xtill'
%     MUD     : fields stored 'Bf','Bm','Bfm','dndt_mud','dBfdt','dBmdt','dBfmdt'
%     i_mc    : index of active coastal element
%
% OUTPUT:
%     COAST     with updated 'mc' (including MUD and DUNE properties)
%     WAVE      with updated 'mc'
%     TRANSP    with updated 'mc'
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

    % temporarily store changes of x,y,xq,yq in a variable to not affect other coastal segments yet.
    fields={'s','ds','xq','yq','PHIc','dPHIc','PHIcxy','PHIcxy0','PHIf','h0','cyclic','clockwise'};
    fieldsdune={'qs','qss','ql','qw','R','SWL','wberm','dfelev','dcelev','xtill'};%,'xhard'};
    fieldsmud={'Bf','Bm','Bfm','dndt_mud','dBfdt','dBmdt','dBfmdt'};
    if DUNE.used
        if isempty(DUNE.xtill)
            fieldsdune=fieldsdune(1:end-1);
        end
        fields={fields{:},fieldsdune{:}};
    end
    if MUD.used
        fields={fields{:},fieldsmud{:}};
    end
    fieldnm{1}=fields;   
    fieldnm{2}={'PHIo','PHItdp','PHIbr','HSo','HStdp','HSbr','dPHIo','dPHItdp','dPHIbr','TP','hbr','dPHIcrit','diff'};
    fieldnm{3}={'QS','QSmax','shadowS','shadowS_h','shadowS_hD','shadow','shadow_h','shadowS_rev','ivals','idrev'}; 
    if i_mc==1
        for kk=1:length(fieldnm{1})
            COAST.([(fieldnm{1}{kk}),'_mc'])=COAST.(fieldnm{1}{kk}); 
            %COAST=rmfield(COAST,fieldnm{1}{kk});
        end
        for kk=1:length(fieldnm{2})
            % add up multiple results of multiple conditions when climate with more than 1 condition at each timestep is evaluated
            if ~isempty(findstr(fieldnm{2}{kk},'PHI')) && kk<=3
                % add up for wave directions
                fldnmHS=regexprep(fieldnm{2}{kk},'PHI','HS');
                HStdp0=max(WAVE.(fldnmHS),1e-6);
                wghtNR=repmat(WAVE.Prob,[1,size(WAVE.(fieldnm{2}{kk}),2)]);
                wghtHS=HStdp0.^2;
                sdr = sum(sind(WAVE.(fieldnm{2}{kk})).*wghtNR.*wghtHS,1);
                cdr = sum(cosd(WAVE.(fieldnm{2}{kk})).*wghtNR.*wghtHS,1);
                WAVE.(fieldnm{2}{kk})=mod(atan2d(sdr,cdr),360);
            elseif size(WAVE.(fieldnm{2}{kk}),1)~=1 
                % add up for wave height of multiple conditions using probability of each wave condition
                wghtNR=repmat(WAVE.Prob,[1,size(WAVE.(fieldnm{2}{kk}),2)]);
                WAVE.(fieldnm{2}{kk})=sum(WAVE.(fieldnm{2}{kk}).*wghtNR,1);
            end
            WAVE.([(fieldnm{2}{kk}),'_mc'])=WAVE.(fieldnm{2}{kk}); 
            %WAVE=rmfield(WAVE,fieldnm{2}{kk});
        end
        for kk=1:length(fieldnm{3})
            if size(TRANSP.(fieldnm{3}{kk}),1)~=1
                % add up for transport of multiple conditions using probability of each wave condition
                TRANSP.(fieldnm{3}{kk})=sum(TRANSP.(fieldnm{3}{kk}).*repmat(WAVE.Prob,[1,size(TRANSP.(fieldnm{3}{kk}),2)]),1);
            end
            TRANSP.([(fieldnm{3}{kk}),'_mc'])=TRANSP.(fieldnm{3}{kk}); 
            %TRANSP=rmfield(TRANSP,fieldnm{3}{kk});
        end
        
    else 
        for kk=1:length(fieldnm{1})
        COAST.([(fieldnm{1}{kk}),'_mc'])=[COAST.([(fieldnm{1}{kk}),'_mc']),nan,COAST.(fieldnm{1}{kk})]; 
        %COAST=rmfield(COAST,fieldnm{1}{kk});
        end
        for kk=1:length(fieldnm{2})
            if size(WAVE.(fieldnm{2}{kk}),1)~=1
                WAVE.(fieldnm{2}{kk})=sum(WAVE.(fieldnm{2}{kk}).*repmat(WAVE.Prob,[1,size(WAVE.(fieldnm{2}{kk}),2)]),1);
            end
            WAVE.([(fieldnm{2}{kk}),'_mc'])=[WAVE.([(fieldnm{2}{kk}),'_mc']),nan,WAVE.(fieldnm{2}{kk})]; 
            %WAVE=rmfield(WAVE,fieldnm{2}{kk});
        end
        for kk=1:length(fieldnm{3})
            if size(TRANSP.(fieldnm{3}{kk}),1)~=1
                TRANSP.(fieldnm{3}{kk})=sum(TRANSP.(fieldnm{3}{kk}).*repmat(WAVE.Prob,[1,size(TRANSP.(fieldnm{3}{kk}),2)]),1);
            end
            TRANSP.([(fieldnm{3}{kk}),'_mc'])=[TRANSP.([(fieldnm{3}{kk}),'_mc']),nan,TRANSP.(fieldnm{3}{kk})]; 
            %TRANSP=rmfield(TRANSP,fieldnm{3}{kk});
        end
    end
    COAST.dSds_mc=[];
    WAVE.diff_mc(isnan(WAVE.diff_mc))=0;
    WAVE.diff_mc=logical(WAVE.diff_mc);
    
    % make sure to store the current grid before merging/splitting in x1_mc, y1_mc, xq1_mc, yq1_mc, PHIcxy1_mc)
    if i_mc==COAST.n_mc && TIME.it==0
        [COAST0]=get_reconnectedgroynes(COAST,GROYNE); 
        COAST.x0_mc=COAST0.x_mc; % x_mc and y_mc have been adjusted in make_sgrid
        COAST.y0_mc=COAST0.y_mc; % x_mc and y_mc have been adjusted in make_sgrid
        COAST.n0_mc=length(find(isnan(COAST.x0_mc)))+1;
    end
    if i_mc==COAST.n_mc 
        COAST.x1_mc=COAST.x_mc; % x_mc and y_mc have been adjusted in make_sgrid
        COAST.y1_mc=COAST.y_mc; % x_mc and y_mc have been adjusted in make_sgrid
        COAST.n1_mc=length(find(isnan(COAST.x1_mc)))+1;
        COAST.xq1_mc=COAST.xq_mc;
        COAST.yq1_mc=COAST.yq_mc;
        COAST.PHIc1_mc=COAST.PHIc_mc;
        COAST.PHIcxy1_mc=COAST.PHIcxy0_mc;
    end
       
    % remove any local fields of the segment that are not further used
    WAVE=rmfield(WAVE,'cbr');
    WAVE=rmfield(WAVE,'nbr');
    WAVE=rmfield(WAVE,'ctdp');
    WAVE=rmfield(WAVE,'ntdp');
    WAVE=rmfield(WAVE,'dPHItdp_sh');
    WAVE=rmfield(WAVE,'dPHIcritbr');
    
    if COAST.i_mc==COAST.n_mc
    WAVE.dPHIcrit_mc0=WAVE.dPHIcrit_mc;
    end
end

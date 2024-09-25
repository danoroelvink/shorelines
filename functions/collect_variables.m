function [COAST,WAVE,TRANSP]=collect_variables(COAST,WAVE,TRANSP,DUNE,MUD,i_mc)
% function [COAST,WAVE,TRANSP]=collect_variables(COAST,WAVE,TRANSP,DUNE,MUD,i_mc)
%
% This routine collects the variables of the evaluation of section 'i_mc'
% into the overarching variable that stores data for all sections, with '_mc' suffix. 
% The x_mc and y_mc are temporarily stored in x1_mc and y1_mc, and written 
% back to x_mc and y_mc when the last element has been evaluated.  
% 
% INPUT:
%     COAST   : fields stored 's','ds','x1','y1','xq','yq','xq1','yq1','PHIc','dPHIc','PHIcxy','PHIf','h0','cyclic','clockwise'
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
    COAST.x1=COAST.x;
    COAST.y1=COAST.y;
    COAST.xq1=COAST.xq;
    COAST.yq1=COAST.yq;
    fields={'s','ds','x1','y1','xq','yq','xq1','yq1','PHIc','dPHIc','PHIcxy','PHIf','h0','cyclic','clockwise'};
    fieldsdune={'qs','qss','ql','qw','R','SWL','wberm','dfelev','dcelev','xtill'};%,'xhard'};
    fieldsmud={'Bf','Bm','Bfm','dndt_mud','dBfdt','dBmdt','dBfmdt'};
    if DUNE.used 
       fields={fields{:},fieldsdune{:}};
    end
    if MUD.used
       fields={fields{:},fieldsmud{:}};
    end
    fieldnm{1}=fields;   
    fieldnm{2}={'PHIo','PHItdp','PHIbr','HSo','HStdp','HSbr','dPHIo','dPHItdp','dPHIbr','TP','hbr','dPHIcrit','diff'};
    fieldnm{3}={'QS','QSmax','shadowS','shadowS_h','shadowS_hD','shadow','shadow_h','ivals','idrev'};
    if i_mc==1
        for kk=1:length(fieldnm{1})
        COAST.([(fieldnm{1}{kk}),'_mc'])=COAST.(fieldnm{1}{kk}); 
        COAST=rmfield(COAST,fieldnm{1}{kk});
        end
        for kk=1:length(fieldnm{2})
        WAVE.([(fieldnm{2}{kk}),'_mc'])=WAVE.(fieldnm{2}{kk}); 
        WAVE=rmfield(WAVE,fieldnm{2}{kk});
        end
        for kk=1:length(fieldnm{3})
        TRANSP.([(fieldnm{3}{kk}),'_mc'])=TRANSP.(fieldnm{3}{kk}); 
        TRANSP=rmfield(TRANSP,fieldnm{3}{kk});
        end
        
    else 
        for kk=1:length(fieldnm{1})
        COAST.([(fieldnm{1}{kk}),'_mc'])=[COAST.([(fieldnm{1}{kk}),'_mc']),nan,COAST.(fieldnm{1}{kk})]; 
        COAST=rmfield(COAST,fieldnm{1}{kk});
        end
        for kk=1:length(fieldnm{2})
        WAVE.([(fieldnm{2}{kk}),'_mc'])=[WAVE.([(fieldnm{2}{kk}),'_mc']),nan,WAVE.(fieldnm{2}{kk})]; 
        WAVE=rmfield(WAVE,fieldnm{2}{kk});
        end
        for kk=1:length(fieldnm{3})
        TRANSP.([(fieldnm{3}{kk}),'_mc'])=[TRANSP.([(fieldnm{3}{kk}),'_mc']),nan,TRANSP.(fieldnm{3}{kk})]; 
        TRANSP=rmfield(TRANSP,fieldnm{3}{kk});
        end
    end
    COAST.dSds_mc=[];
    WAVE.diff_mc(isnan(WAVE.diff_mc))=0;
    WAVE.diff_mc=logical(WAVE.diff_mc);
    
    % make sure to get rid of the temporary x1,y1,xq1,yq1
    % and put the x1_mc,y1_mc,xq1_mc,yq1_mc back to the regular mc
    if i_mc==COAST.n_mc
        COAST.x_mc=COAST.x1_mc;
        COAST.y_mc=COAST.y1_mc;
        COAST.xq_mc=COAST.xq1_mc;
        COAST.yq_mc=COAST.yq1_mc;
    end
       
    % remove any local fields of the segment that are not further used
    WAVE=rmfield(WAVE,'cbr');
    WAVE=rmfield(WAVE,'nbr');
    WAVE=rmfield(WAVE,'ctdp');
    WAVE=rmfield(WAVE,'ntdp');
    WAVE=rmfield(WAVE,'dPHItdp_sh');
    WAVE=rmfield(WAVE,'dPHIcritbr');
    WAVE=rmfield(WAVE,'diff2');
    
    if COAST.i_mc==COAST.n_mc
    WAVE.dPHIcrit_mc0=WAVE.dPHIcrit_mc;
    end
end

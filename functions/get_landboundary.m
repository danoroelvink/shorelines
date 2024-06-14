function [LDB]=get_landboundary(LDBfilename,closepoly)
%function [LDB]=get_landboundary(LDBfilename,closepoly)
%
% function read LDBs and POLs
% also more than one polygon can be read.
%
% INPUT:
%    LDBfilename   String with ldb filename
%    closepoly     (optional) Switch (0/1) for either making sure to close each polygon or not (by means of the last point being a copy of the first)
%                      closepoly=0 (default when not specified) : using one on one only the data that are in the file (irrespectyive of whether it is open or not)
%                      closepoly=-1 : make open line sections (i.e. remove the last point of each polygon if it is copy of the first)
%                      closepoly=1  : make closed polygons (i.e. add the last point of each polygon as a copy of the first point, if it is not already the case)
% 
% OUTPUT:
%    LDB           Coordinates of the landboundary file
%
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2022 IHE Delft & Deltares
%
%       Bas Huisman
%       bas.huisman@deltares.nl
%       Boussinesqweg 1
%       2629HV Delft
%
%       Dano Roelvink
%       d.roelvink@un-ihe.org
%       Westvest 7
%       2611AX Delft
%
%   This library is free software: you can redistribute TIME.it and/or
%   modify TIME.it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation, either
%   version 2.1 of the License, or (at your option) any later version.
%
%   This library is distributed in the hope that TIME.it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses
%   --------------------------------------------------------------------

    LDB          = [];
    defaultNAN   = 999.999;
    jj           = 0;
    linenr       = 0;
    if nargin==1
    closepoly    = 0;
    end
    
    %% OPEN FILE
    fid = fopen(LDBfilename,'r');
    if fid==-1
        fprintf('Error : ''%s'' does not exist!\n',LDBfilename)
    else
        %% READ FILE
        while ~feof(fid)
            line = fgetl(fid); line=[line,' ']; 
            if isempty(findstr(line(1),'*'));
                %% READ HEADER
                linehdr = line;
                line = fgetl(fid); 
                if ~isempty(deblank(line)) && ~feof(fid)
                    [nrpoints,nrcolumns]=strread([line,'  '],'%f %f','delimiter',' ');
                    if jj>0
                        jj=jj+1;
                        LDB(jj,1) = nan;
                        LDB(jj,2) = nan;
                        if nrcolumns>2
                            LDB(jj,3) = nan;
                        end
                    end
                    
                    %% READ DATA
                    jj0=jj;
                    for ii=1:nrpoints
                        line = fgetl(fid); 
                        [a,b,c]=strread([line,'  '],'%f %f %f','delimiter',' ');
                        jj=jj+1;
                        LDB(jj,1) = a;
                        LDB(jj,2) = b;
                        if nrcolumns>2
                            LDB(jj,3) = c;
                        end
                    end
                        
                    %% OPEN POLYGONS (if -1)
                    if closepoly==-1 && jj>jj0+1
                        if LDB(jj,1)==LDB(jj0+1,1) && LDB(jj,2)==LDB(jj0+1,2)
                            LDB=LDB(1:jj-1,:);
                            jj=jj-1;
                        end
                    %% CLOSE POLYGONS (if 1)
                    elseif closepoly==1 && jj>jj0+1
                        if LDB(jj,1)~=LDB(jj0+1,1) || LDB(jj,2)~=LDB(jj0+1,2)
                            LDB=[LDB;LDB(jj0+1,:)];
                            jj=jj+1;
                        end
                    end            
                        
                    % replace NAN's
                    LDB(LDB==defaultNAN)=nan;
                end
            end
        end
        fclose(fid);
        
    end
end

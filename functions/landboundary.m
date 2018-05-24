function varargout=landboundary(cmd,varargin)
%LANDBOUNDARY Read/write land boundary files.
%   XY = LANDBOUNDARY('read',FILENAME) reads the specified file and returns
%   the data as one Nx2 array.
%
%   [X,Y] = LANDBOUNDARY(...) returns separate X and Y arrays.
%
%   LANDBOUNDARY('write',FILENAME,XY) writes a landboundary to file. XY
%   should either be a Nx2 array containing NaN separated line segments or
%   a cell array containing one line segment per cell.
%
%   LANDBOUNDARY('write',FILENAME,X,Y) writes a landboundary to file. X
%   and Y supplied as separate Nx1 arrays containing NaN separated line
%   segments or cell arrays containing one line segment per cell. The X and
%   Y line segments should correspond in length.
%
%   LANDBOUNDARY(...,'-1') does not write line segments of length 1.
%
%   LANDBOUNDARY(...,'dosplit') saves line segments as separate TEKAL
%   blocks instead of saving them as one long line interrupted by missing
%   values. This approach is well suited for spline files, but less suited
%   for landboundaries with a large number of segments.
%
%   FILE = LANDBOUNDARY('write',...) returns a file info structure for the
%   file written. This structure can be used to read the file using the
%   TEKAL function.
%
%   See also TEKAL.

%----- LGPL --------------------------------------------------------------------
%                                                                               
%   Copyright (C) 2011-2012 Stichting Deltares.                                     
%                                                                               
%   This library is free software; you can redistribute it and/or                
%   modify it under the terms of the GNU Lesser General Public                   
%   License as published by the Free Software Foundation version 2.1.                         
%                                                                               
%   This library is distributed in the hope that it will be useful,              
%   but WITHOUT ANY WARRANTY; without even the implied warranty of               
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
%   Lesser General Public License for more details.                              
%                                                                               
%   You should have received a copy of the GNU Lesser General Public             
%   License along with this library; if not, see <http://www.gnu.org/licenses/>. 
%                                                                               
%   contact: delft3d.support@deltares.nl                                         
%   Stichting Deltares                                                           
%   P.O. Box 177                                                                 
%   2600 MH Delft, The Netherlands                                               
%                                                                               
%   All indications and logos of, and references to, "Delft3D" and "Deltares"    
%   are registered trademarks of Stichting Deltares, and remain the property of  
%   Stichting Deltares. All rights reserved.                                     
%                                                                               
%-------------------------------------------------------------------------------
%   http://www.deltaressystems.com
%   $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/trunk/src/tools_lgpl/matlab/quickplot/progsrc/landboundary.m $
%   $Id: landboundary.m 1147 2011-12-31 23:43:35Z jagers $

if nargout>0
    varargout=cell(1,nargout);
end
if nargin==0
    return
end
switch cmd
    case 'read'
        Out=Local_read_file(varargin{:});
        if nargout==1
            varargout{1}=Out;
        elseif nargout>1
            varargout{1}=Out(:,1);
            varargout{2}=Out(:,2);
        end
    case 'write'
        Out=Local_write_file(varargin{:});
        if nargout>0
            varargout{1}=Out;
        end
    otherwise
        uiwait(msgbox('unknown command','modal'));
end


function Data=Local_read_file(filename)
Data=[];
if nargin==0
    [fn,fp]=uigetfile('*.ldb');
    if ~ischar(fn)
        return
    end
    filename=[fp fn];
end

T=tekal('open',filename);

lasterr('')
try
    Sz=cat(1,T.Field.Size);
    if ~all(Sz(:,2)==Sz(1,2))
        error('The number of columns in the files is not constant.')
    end
    Sz=[sum(Sz(:,1))+size(Sz,1)-1 Sz(1,2)];
    offset=0;
    Data=repmat(NaN,Sz);
    for i=1:length(T.Field)
        Data(offset+(1:T.Field(i).Size(1)),:)=tekal('read',T,i);
        offset=offset+T.Field(i).Size(1)+1;
    end
    Data( (Data(:,1)==999.999) & (Data(:,2)==999.999) ,:)=NaN;
catch
    fprintf(1,'ERROR: Error extracting landboundary from tekal file:\n%s\n',lasterr);
end


function TklFileInfo=Local_write_file(filename,varargin)

if nargin==1
    [fn,fp]=uiputfile('*.*');
    if ~ischar(fn)
        TklFileInfo=[];
        return
    end
    filename=[fp fn];
end

j=0;
RemoveLengthOne=0;
XYSep=0;
DoSplit=0;
CellData=0;
Format='%4i';
i=1;
while i<=nargin-1
    if ischar(varargin{i}) && strcmp(varargin{i},'-1')
        RemoveLengthOne=1;
    elseif ischar(varargin{i}) && strcmpi(varargin{i},'dosplit')
        DoSplit=1;
    elseif ischar(varargin{i}) && strcmpi(varargin{i},'format') && i<nargin-1
        Format=varargin{i+1};
        i=i+1;
    elseif isnumeric(varargin{i}) && (j==0 || j==1)
        if j==0
           data1=varargin{i};
           j=1;
        else % j==1
            data2=varargin{i};
            XYSep=1; % x and y supplied separately?
        end
    elseif iscell(varargin{i}) && (j==0 || j==1)
        CellData=1;
        if j==0
           Data1=varargin{i};
           j=1;
        else % j==1
            Data2=varargin{i};
            XYSep=1; % x and y supplied separately?
        end
    else
        error('Invalid input argument %i',i+2)
    end
    i=i+1;
end

if CellData
    Ncell = numel(Data1);
else
    Ncell = 1;
end

j=0;
for c = 1:Ncell
    if CellData
        data1 = Data1{c};
        if XYSep
            data2 = Data2{c};
            if size(data1,1)~=numel(data1)
                data1 = data1(:);
                data2 = data2(:);
            end
        else
            if ndims(data1)>2
                error('Invalid size of data array %i.',c)
            elseif size(data1,2)~=2
                if size(data1,1)==2
                    data1 = data1';
                else
                    error('Invalid size of data array %i.',c)
                end
            end
        end
    end
    %
    I=[0; find(isnan(data1(:,1))); size(data1,1)+1];
    if ~DoSplit
        I=[0;size(data1,1)+1];
        data1(isnan(data1))=999.999;
        if XYSep
            data2(isnan(data2))=999.999;
        end
    end
    for i=1:(length(I)-1)
        if I(i+1)>(I(i)+1+RemoveLengthOne)
            % remove lines of length 0 (and optionally 1)
            j=j+1;
            if j==1
                T.Field(1).Comments = { ...
                    '*column 1 = x coordinate'
                    '*column 2 = y coordinate'};
            end
            T.Field(j).Name = sprintf(Format,j);
            T.Field(j).Size = [I(i+1)-I(i)-1 2];
            if XYSep
                T.Field(j).Data = [data1((I(i)+1):(I(i+1)-1)) data2((I(i)+1):(I(i+1)-1))];
            else
                T.Field(j).Data = data1((I(i)+1):(I(i+1)-1),:);
            end
        end
    end
end

TklFileInfo=tekal('write',filename,T);

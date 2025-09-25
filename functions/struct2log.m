function struct2log(P,structname,writeopt)
% function struct2log(P,structname,writeopt)
% 
% This function logs the input parameters of a model run to 'logfile.txt'
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2024 IHE Delft & Deltares
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

    % Open a file for appending
    fid = fopen('logfile.txt', writeopt);

    if ischar(P)
        % print the characters of a warning message to the log-file
        if ischar(structname)
            structname={structname};
        end
        for jj=1:length(structname)
            fprintf(fid, '%s\n',[structname{jj}]);
            if strcmpi(P,'print') || strcmpi(P,'warning')
                fprintf('%s\n',[structname{jj}]);
            end
        end
    else
        % Get the field names of the P structure
        fieldname = fieldnames(P);
        % Write the structure contents to the file
        fprintf(fid, '%s\n',['%%%%%%%%%%% ',structname,' ',repmat('%',[1,15-length(structname)])]);
        for i=1:length(fieldname)
            val = P.(fieldname{i});
            fieldnm = fieldname{i};
            printinputfield(fid,fieldnm,val,' ');
        end
    end
    fprintf(fid, '\n');
    % Close the file
    fclose(fid);

end

function printinputfield(fid,fieldnm,val,addstr,endline)

    endofline='\n';
    nrspaces=26;
    if nargin==5
        endofline='';
        nrspaces=1;
        fieldnm=['',fieldnm,repmat(' ',[1 max(nrspaces-length(fieldnm),0)])];
    else
        fieldnm=['',fieldnm,repmat(' ',[1 max(nrspaces-length(fieldnm),0)]),' ='];
    end
    % loop over content

    % fill out fieldname with blank spaces
    

    % empty field
    if isempty(val)
       fprintf(fid, ['%s%s %s',endofline],addstr,fieldnm,'[]');
  
    % numeric / vector / matrix input
    elseif isnumeric(val)
        if length(val)==1
            fprintf(fid, ['%s%s %s',endofline],addstr,fieldnm,num2str(val(:)'));
        else
            sz=size(val);
            if sz(2)==1 && sz(1)>1
               val=val';
               sz=size(val);
            end
            fprintf(fid, '%s%s [',addstr,fieldnm);
            for mm=1:sz(1)
               for nn=1:sz(2)               
               fprintf(fid, ' %s', num2str(val(mm,nn)'));
               end
               if mm<sz(1)
               fprintf(fid, [' ... \n',repmat(' ',[1 nrspaces+5])]);
               end
            end
            fprintf(fid, [' ] ',endofline]);
        end
        
    % cell input
    elseif iscell(val)
       sz=size(val);
       if sz(2)==1 && sz(1)>1
           val=val';
           sz=size(val);
       end
       fprintf(fid, '%s%s {',addstr,fieldnm);
       for mm=1:sz(1)
           for nn=1:sz(2)               
           %fprintf(fid, ' %s', num2str(val{mm,nn}'));
           fieldnm2='';
           printinputfield(fid,fieldnm2,val{mm,nn}','',1);
           end
           if mm<sz(1)
           fprintf(fid, [' ... \n',repmat(' ',[1 nrspaces+5])]);
           end
       end
       fprintf(fid, [' } ',endofline]);

    % structures
    elseif isstruct(val)
        fieldnm=deblank(fieldnm);
        fields2=fields(val);
        %fprintf(fid, '%s%s   %s\n',addstr,fieldnm,' % <struct>');
        for mm=1:length(val)
        for gg=1:length(fields2)
           %fprintf(fid, '%32s%s\n','', ['.',fields2{gg}]);
           addstr=' ';
           if length(val)==1
               fieldnm2=['',[deblank(fieldnm(1:end-1)),'.',fields2{gg}],repmat(' ',[1 max(nrspaces-length([deblank(fieldnm(1:end-1)),'.',fields2{gg}]),0)])];
           else
               fieldnm2=['',[deblank(fieldnm(1:end-1)),'(',num2str(mm),').',fields2{gg}],repmat(' ',[1 max(nrspaces-length([deblank(fieldnm(1:end-1)),'(',num2str(mm),').',fields2{gg}]),0)])];
           end
           printinputfield(fid,fieldnm2,val(mm).(fields2{gg}),addstr);
        end
        end
       
    % characters / text
    elseif ischar(val)
        if isempty(addstr)
            fprintf(fid, ['%s%s%s',endofline],addstr,fieldnm,['''',val(:)','''']);
        else
            try
              fprintf(fid, ['%s%s %s',endofline],addstr,fieldnm,['''',val(:)','''']);
            catch
              fprintf(fid, ['%s%s ',endofline],addstr,['''',fieldnm,'''']);
            end
        end
    end
end
function struct2log(P,structname,writeopt)

    % Get the field names of the P structure
    fieldname = fieldnames(P);
    % Open a file for appending
    fid = fopen('logfile.txt', writeopt);

    % Write the structure contents to the file
    fprintf(fid, '%s\n',['%%%%%%%%%%% ',structname,' ',repmat('%',[1,15-length(structname)])]);
    for i=1:length(fieldname)
        val = P.(fieldname{i});
        fieldnm = fieldname{i};
        printinputfield(fid,fieldnm,val,' ');
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
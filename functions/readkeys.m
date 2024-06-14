function [data] = readkeys(filename)
% function [data] = readkeys(filename)
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

   number=0;
   fid=fopen(filename);
   while 1
      s=fgetl(fid);
      if ~isempty(s)
         if ~ischar(s), break, end
         s = strip(s);
         cont=strfind(s,'...');
         while ~isempty(cont)
            nxt=fgetl(fid);
            s=[s(1:cont-1),strip(nxt)];
            cont=strfind(s,'...');
         end
         if s(1)~='%'
            k=find(s=='=',1);
            if ~isempty(k)
               number=number+1;
               keywords{number}=strip(s(1:k-1));
               val=strip(s(k+1:end));
               % m=find(val==' '|val==';'|val=='%'|val=='/',1);
               %$m=find(val==';'|val=='%'|val=='/',1);
               m=find(val=='%'|val=='/',1);
               if isempty(m)
                  values{number}=val;
               else
                  values{number}=val(1:m-1);
               end
            end
         end
      end
   end
   fclose(fid);
   for i=1:length(keywords);
      keyword=keywords{i};
      value=values{i};
      eval(['data.',keyword,'=',value,';']);
   end
end

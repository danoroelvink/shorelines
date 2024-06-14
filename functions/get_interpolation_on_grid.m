function [var1i,var2i,idGRID,distw]=get_interpolation_on_grid(method,xq,yq,xw,yw,var1,var2,fldnms1,fldnms2)
% function [var1i,var2i,idGRID,distw]=get_interpolation_on_grid(method,xq,yq,xw,yw,var1,var2,fldnms1,fldnms2);
% function [var1i,idGRID,distw]=get_interpolation_on_grid(method,xq,yq,xw,yw,var1,fldnms1,fldnms2);
%
% find the right alongshore location for each of the wave climates (xq and yq = transport points)
% sort locations alongshore
% choose only the closest wave climate at each grid cell, to make sure that no multiple climates are enforced on a single grid cell.
% Make sure that wave climates are not too far from the grid cells. 
% Otherwise throw out the ones that are at great distance.
%  - remove wave climate points that are further away (in distance) 
%  - from the grid point than the distance to another nearby climate point.
%  - and which are more than 2x the cross-shore distance from the coast than the adjacent climate point
% interpolate the waves at the right alongshore location 
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2020 IHE Delft & Deltares
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

    % allow model to deal with situation where 'var1' and/or 'var2' are specified as a vector. 
    % and the situation wherein 'var1' and/or 'var2' are not specified. 
    if nargin<6
        var1=struct;
    elseif isempty(var1)
        var1=struct;
    elseif isnumeric(var1)
        var1tmp=var1;
        var1=struct;
        var1.data=var1tmp;
    elseif isstruct(var1) && nargin==8
        fldnms2=fldnms1;
        fldnms1=var2;
        %var1=struct;
        %for ff=1:length(fldnms2)
        %    var1.(fldnms2{ff})=var1.(fldnms2{ff});
        %end
    elseif isstruct(var1) && nargin==9
        var1tmp=var1;
        var1=struct;
        for ff=1:length(fldnms1)
            var1.(fldnms1{ff})=var1tmp.(fldnms1{ff});
        end
    elseif isstruct(var1)
        % default, use all fields
    end
    if nargin<7
        var2=struct;      
    elseif isempty(var2)
        var2=struct;
    elseif isnumeric(var2)
        var2tmp=var2;
        var2=struct;
        var2.data=var2tmp;
    elseif isstruct(var1) && nargin==8
        var2=struct;
        for ff=1:length(fldnms2)
            var2.(fldnms2{ff})=var1.(fldnms2{ff});
        end
    elseif isstruct(var2) && nargin==9
        % create a new structure with the specified fields specified structure with the new interpolated data
        var2tmp=var2;
        var2=struct;
        for ff=1:length(fldnms2)
            var2.(fldnms2{ff})=var2tmp.(fldnms2{ff});
        end
    elseif isstruct(var2)
        % default, use all fields
    end

    % find the fields with the data
    % the scalars as fields in 'var1'
    % any field with the directions in degrees in 'var2'
    if nargin<8
    fldnms1=get_fields(var1);
    fldnms2=get_fields(var2);
    end
    ff1=length(fldnms1);
    ff2=length(fldnms2);
    var1i=struct;
    var2i=struct;
    idGRID=[];      % these are either the grid locations for each of the applied wave climates 'WVC(kk)' in case of 'alongshore_mapping' or the applied wave cliamte points in sorted-order for each grid cell 'weighted_distance' (with a NaN meaning the wave climate file could not be used, e.g. too far from coastline)
    distw=[];       % distance used for the shadowing function for the length of the ray (for the offshore part)
    
    % dealing with larger data fields [MxN] instead of [1xN]
    transposedata=zeros(ff1+ff2,1);
    datalength=ones(ff1+ff2,1);
    for kk=1:length(fldnms1)
        if size(var1.(fldnms1{kk}),2)~=length(xw) %&& size(var1.(fldnms1{kk}),2)~=1
            var1.(fldnms1{kk})=var1.(fldnms1{kk})';
            transposedata(kk)=1;
        end
        if size(var1.(fldnms1{kk}),1)~=length(xw) %&& size(var1.(fldnms1{kk}),1)~=1
            nm=size(var1.(fldnms1{kk}),1);
            datalength(kk)=nm;
        end
    end
    for kk2=1:length(fldnms2)
        if size(var2.(fldnms2{kk2}),2)~=length(xw) %&& size(var2.(fldnms2{kk2}),2)~=1
            var2.(fldnms2{kk2})=var2.(fldnms2{kk2})';
            transposedata(kk+kk2)=1;
        end
        if size(var2.(fldnms2{kk2}),1)~=length(xw) %&& size(var2.(fldnms2{kk2}),1)~=1
            nm=size(var2.(fldnms2{kk2}),1);
            datalength(kk+kk2)=nm;
        end
    end
    
    % catch input of method if it is defined slightly differently
    if ~isempty(findstr(method,'weight')) || ~isempty(findstr(method,'distance'))
        method='weighted_distance';
    elseif ~isempty(findstr(method,'alongshore')) || ~isempty(findstr(method,'mapping'))
        method='alongshore_mapping';
    end
    
    %------------------------------------------------------------------------%
    %% use quadratic distance weighting for each of the wave locations
    %------------------------------------------------------------------------%          
    if strcmpi(method,'weighted_distance') && length(xw)>=1

        pwr=1;        % use linear distance weighting
        nrpoints=2;   % number of nearest points used (rest at larger distance is omitted)
        for ii=1:length(xq)
            dist=((xq(ii)-xw).^2+(yq(ii)-yw).^2).^0.5;
            dist=dist(:)';
            [dists,ids]=sort(dist);
            
            if length(ids)>nrpoints
                dists=dists(1:nrpoints);
                ids=ids(1:nrpoints);
            end
            distw(1,ii)=min(dists); % distance used for the shadowing function for the length of the ray

            if length(ids)>=2
                wght0=1-(dists/max(sum(dists),eps));
                wght=wght0.^pwr/max(sum(wght0.^pwr),eps);
                %fprintf('ii=%3.0f : ',ii);
                %fprintf('dists=%4.3f, ',dists);
                %fprintf('ids=%4.3f, ',ids);
                %fprintf('wghts=%4.3f, ',wght);
                %fprintf('\n');
                for kk=1:length(fldnms1)
                    for gg=1:datalength(kk)
                    var1i.(fldnms1{kk})(gg,ii)=sum(var1.(fldnms1{kk})(gg,ids).*wght);
                    end
                end
                for kk=1:length(fldnms2)
                    for gg=1:datalength(ff1+kk)
                    cphi=sum(cosd(var2.(fldnms2{kk})(gg,ids)).*wght);
                    sphi=sum(sind(var2.(fldnms2{kk})(gg,ids)).*wght);
                    var2i.(fldnms2{kk})(gg,ii)=mod(atan2d(sphi,cphi),360);
                    end
                end
                idGRID(ii,1)=[(1-wght(1,1))*diff(ids(1:2))+ids(1)]; % ids(:)';
            elseif length(ids)==1
                for kk=1:length(fldnms1)
                    for gg=1:datalength(kk)
                    var1i.(fldnms1{kk})(gg,ii)=var1.(fldnms1{kk})(gg,ids);
                    end
                end
                for kk=1:length(fldnms2)
                    for gg=1:datalength(ff1+kk)
                    var2i.(fldnms2{kk})(gg,ii)=mod(var2.(fldnms2{kk})(gg,ids),360);
                    end
                end
                idGRID(ii,1)=ids; 
            else
                for kk=1:length(fldnms1)
                    var1i.(fldnms1{kk})=zeros(datalength(kk),length(xq));
                end
                for kk=1:length(fldnms2)
                    var2i.(fldnms2{kk})=zeros(datalength(ff1+kk),length(xq));
                end
            end
        end
    
    %------------------------------------------------------------------------%
    %% Map the wave climate locations on the coast and interpolate alongshore 
    %------------------------------------------------------------------------%
    elseif strcmpi(method,'alongshore_mapping') && length(xw)>=1
        
        % find the right alongshore location for each of the wave climates
        idGRID=[];
        distq = [0,cumsum((diff(xq).^2+diff(yq).^2).^0.5)];
        dist0=[];
        distw=nan([1,length(xq)]);
        for kk=1:length(xw)
            dist=((xq-xw(kk)).^2+(yq-yw(kk)).^2).^0.5;
            id=find(dist==min(dist),1);
            distxy=distq(id);
            idGRID(kk,1)=id;            
            dist0(kk,1)=distxy;
            distw(1,id)=min(dist); 
        end
        distw=interpNANs(distw); % distance used for the shadowing function for the length of the ray (for the offshore part)
        idGRID0=idGRID;
        
        % sort locations alongshore
        [~,idu]=unique(dist0);
        [WVCsort,IDS]=sort(idGRID);
        IDS=intersect(IDS,idu);       
        dist0=dist0(IDS);
        idGRID=idGRID(IDS); 
        for kk=1:length(fldnms1)
            var1.(fldnms1{kk})=var1.(fldnms1{kk})(:,IDS(:)');
        end
        for kk=1:length(fldnms2)
            var2.(fldnms2{kk})=var2.(fldnms2{kk})(:,IDS(:)');
        end
        xw=xw(IDS);
        yw=yw(IDS);
        
        % choose only the closest wave climate at each grid cell, 
        % to make sure that no multiple climates are enforced on a single grid cell.
        dist3=[];
        idDATA=[1:length(idGRID)];
        ids2=[];
        for xx=1:length(idGRID)
            ID=find(idGRID==idGRID(xx));
            dist2=[];
            for dd=1:length(ID)
                dist1=((xq-xw(ID(dd))).^2+(yq-yw(ID(dd))).^2).^0.5;
                dist2(dd)=min(dist1);
            end
            idDATA=setdiff(idDATA,ID(find(dist2~=min(dist2))));
            ids2=[ids2;ID(find(dist2~=min(dist2)))];
            dist3(xx)=min(dist2);
        end
        dist0=dist0(idDATA);
        for kk=1:length(fldnms1)
            var1.(fldnms1{kk})=var1.(fldnms1{kk})(:,idDATA(:)');
        end
        for kk=1:length(fldnms2)
            var2.(fldnms2{kk})=var2.(fldnms2{kk})(:,idDATA(:)');
        end
        idGRID=idGRID(idDATA);
        idGRID0(IDS(ids2))=nan;
        
        % Make sure that wave climates are not too far from the grid cells. 
        % Otherwise throw out the ones that are at great distance.
        iter=10;
        IDS0=IDS;
        ids3=[];
        for kk=1:iter 
            dx=[0,cumsum((diff(xq).^2+diff(yq).^2).^0.5)];
            dx0=((xq(end)-xq(1)).^2+(yq(end)-yq(1)).^2).^0.5;
            dxend=dx(end)-dist0(end);
            dx1=dist0(1)+dxend;
            dx=[dx0+dx1,diff(dist0(:))',dx0+dx1];                          % alongshore distance between grid points of 2 consecutive climates
            dcross=[dist3(1)-dist3(end),diff(dist3),dist3(1)-dist3(end)];  % difference in nearest cross-shore distance to wave climates of 2 consecutive climates
            dist3=[dist3(end),dist3,dist3(1)];
            idxx=[1:length(idGRID)];
            for xx=1:length(idGRID)
                % remove wave climate points that are further away (in distance) 
                % from the grid point than the distance to another nearby climate point.
                % and which are more than 2x the cross-shore distance from the coast than the adjacent climate point
                distfactor=2;
                if (dcross(xx)>dx(xx) && dist3(xx+1)>distfactor*dist3(xx)) || ...     % forward check if cross-shore distance to wave output location is not increasing excessively
                   (-dcross(xx+1)>dx(xx+1) && dist3(xx+1)>distfactor*dist3(xx+2))     % backward check if cross-shore distance to wave output location is not increasing excessively
                    idxx=setdiff(idxx,xx);
                    ids3=[ids3;IDS0(xx)];
                end
            end
            dist0=dist0(idxx);
            dist3=dist3(idxx+1);
            idGRID=idGRID(idxx);
            for kk=1:length(fldnms1)
                for gg=1:datalength(kk)
                var1.(fldnms1{kk}(gg,:))=var1.(fldnms1{kk})(gg,idxx(:)');
                end
            end
            for kk=1:length(fldnms2)
                for gg=1:datalength(ff1+kk)
                var2.(fldnms2{kk}(gg,:))=var2.(fldnms2{kk})(gg,idxx(:)');
                end
            end
            IDS0=IDS0(idxx);
        end
        
        % interpolate the waves at the right alongshore location 
        if length(dist0)>=2
            for kk=1:length(fldnms1)
                for gg=1:datalength(kk)
                var1i.(fldnms1{kk})(gg,:)=interp1(dist0,var1.(fldnms1{kk})(gg,:),distq,'linear');
                var1i.(fldnms1{kk})(gg,:)=interpNANs(var1i.(fldnms1{kk})(gg,:));
                end
            end
            for kk=1:length(fldnms2)
                for gg=1:datalength(ff1+kk)
                cphi=interp1(dist0,cosd(var2.(fldnms2{kk})(gg,:)),distq,'linear');
                sphi=interp1(dist0,sind(var2.(fldnms2{kk})(gg,:)),distq,'linear');
                cphi=interpNANs(cphi);
                sphi=interpNANs(sphi);
                var2i.(fldnms2{kk})(gg,:)=mod(atan2d(sphi,cphi),360);
                end
            end
            idGRID(ids3)=nan; 
        elseif length(dist0)==1
            for kk=1:length(fldnms1)
                for gg=1:datalength(kk)
                var1i.(fldnms1{kk})(gg,:)=repmat(var1.(fldnms1{kk})(gg,:),[1,length(distq)]);
                end
            end
            for kk=1:length(fldnms2)
                for gg=1:datalength(ff1+kk)
                var2i.(fldnms2{kk})(gg,:)=repmat(mod(var2.(fldnms2{kk})(gg,:),360),[1,length(distq)]);
                end
            end
            idGRID(ids3)=nan; 
        else
            for kk=1:length(fldnms1)
                var1i.(fldnms1{kk})=zeros(datalength(kk),length(xq));
            end
            for kk=1:length(fldnms2)
                var2i.(fldnms2{kk})=zeros(datalength(ff1+kk),length(xq));
            end
        end
    end
    
    % transpose data back
    for kk=1:length(fldnms1)
        if transposedata(kk)==1
            var1i.(fldnms1{kk})=var1i.(fldnms1{kk})';
        end
    end
    for kk=1:length(fldnms2)
        if transposedata(ff1+kk)==1
            var2i.(fldnms2{kk})=var2i.(fldnms2{kk})';
        end
    end

    % re-organize output data back to a vector if it was specified as a vector in the input
    if exist('var1tmp') && isfield(var1i,'data')
        var1i=var1i.data;
    end
    if exist('var2tmp') && isfield(var2i,'data')
        var2i=var2i.data;
    end    
    
    % in case the input is only for var1 and not for var2, then combine both.
    if nargin==8
        for kk=1:length(fldnms2)
            var1i.(fldnms2{kk})=var2i.(fldnms2{kk});
        end
        var2i=idGRID;
        idGRID=distw;
    end
end

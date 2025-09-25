function [NOUR]=get_nourishments(TIME,COAST,NOUR)
% function [NOUR]=get_nourishments(TIME,COAST,NOUR)
%
% INPUT: 
%    TIME
%       .tnow            : current date in the model (in timenum format)
%    COAST
%       .x               : x-coordinates of the coastline [m]
%       .y               : y-coordinates of the coastline [m]
%       .n               : Number of coastline points (length of x and y)
%    NOUR
%       .nourish         : switch for using nourishments (0, 1 or 2) -> 2 is common method used
%       .growth	         : factor for scaling the nourishment efficiency (from 0 to 1, default is 1)
%       .xnour           : x-coordinates of nourishments [Nx1]
%       .ynour           : y-coordinates of nourishments [Nx1]
%       .nnour           : number of nourishments
%       .tstart          : start time of nourishments [days in datenum format]
%       .tend            : end time of nourishments [days in datenum format]
%       .rate            : nourishment rate [m3/yr]
%
% OUTPUT:
%     NOUR                    
%       .rate_density    : actual nourishment rate in [m3/m/year]
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

    NOUR.rate_density=zeros(1,COAST.n);
    if NOUR.nourish==1
        for i_n=1:NOUR.nnour
            if TIME.tnow>=NOUR.tstart(i_n) && TIME.tnow<=NOUR.tend(i_n)
                [ x_n,y_n,~,~,~ ] = get_one_polygon( NOUR.xnour,NOUR.ynour,i_n );
                in=inpolygon(COAST.x,COAST.y,x_n,y_n);
                if sum(in)>0
                   cl=(hypot(diff(COAST.x(in)),diff(COAST.y(in))));
                   idnotnan=find(~isnan(cl));
                   clsum=sum(cl(idnotnan));
                   NOUR.rate_density=NOUR.rate_density+in*NOUR.rate(i_n)/clsum;
                end
            end
        end
    elseif NOUR.nourish==2
        for i_n=1:NOUR.nnour
            if TIME.tnow>=NOUR.tstart(i_n) && TIME.tnow<=NOUR.tend(i_n)
                if ~isempty(findstr(lower(NOUR.method),'default'))

                    %% default method : just use begin and end point of nourishments
                    % find indices of begin and end point of each nourishment
                    % compute number of grid cells in-between begin and end point
                    dxnour=[];
                    iddist=[];
                    idx={};
                    for i_mc=1:COAST.n_mc
                        [ xc,yc,~,~,~ ] = get_one_polygon( COAST.x_mc,COAST.y_mc,i_mc );
                        dist1=((xc-NOUR.xnour(i_n,1)).^2+(yc-NOUR.ynour(i_n,1)).^2).^0.5;
                        dist2=((xc-NOUR.xnour(i_n,2)).^2+(yc-NOUR.ynour(i_n,2)).^2).^0.5;
                        dist3=((xc-mean(NOUR.xnour(i_n,:))).^2+(yc-mean(NOUR.ynour(i_n,:))).^2).^0.5;
                        idx1=find(dist1==min(dist1),1);
                        idx2=find(dist2==min(dist2),1);
                        idx3=find(dist3==min(dist3),1);
                        idx{i_mc}=[min([idx1,idx2,idx3]):max([idx1,idx2,idx3])]; % find the coastline id's that are affected
                        iddist(i_mc,:)=[min(dist1),min(dist2),min(dist3)]; % store the minimum distance to this coastal section
                    
                        % Compute grid cell size
                        s=[0,cumsum((diff(xc).^2+diff(yc).^2).^0.5)];
                        ds0=diff(s);
                        ds=[ds0(1),(ds0(1:end-1)+ds0(2:end))/2,ds0(end)]; 
                        % Check if coastal section is cyclic
                        cyclic=get_cyclic(xc,yc,COAST.ds0);
                        if cyclic
                            ds(1)=(ds0(1)+ds0(end))/2;
                            ds(end)=(ds0(1)+ds0(end))/2; 
                        end
                        % determine the total length of coastline that is nourished per coastal element
                        dxnour(i_mc)=sum(ds(idx{i_mc})); 
                        if isnan(dxnour(i_mc))
                            dxnour(i_mc)=0;
                        end
                    end 
                    % don't nourish coastal sections where there is only 1 nearby point
                    % and which is further away from nourishment than another coastal section
                    % also remove elements where not any of the end-points is closest to that particular coastal elements (but closer to another coastline element)
                    for i_mc=1:COAST.n_mc
                        if length(idx{i_mc})>=1 && iddist(i_mc,1)>min(iddist(:,1)) && iddist(i_mc,2)>min(iddist(:,2)) && iddist(i_mc,3)>min(iddist(:,3))
                            idx{i_mc}=[];
                            dxnour(i_mc)=0;
                        end
                    end
                else

                    %% Alternative method : chop the nourishment up in 20 points along its length
                    % this method can be used by setting S.nourmethod to soemthing else than 'default'
                    % it is a slower method, but is better suitable when multiple coastline elements are used (e.g. with multiple islands)
                    dist0=[0,(diff(NOUR.xnour(i_n,:)).^2+diff(NOUR.ynour(i_n,:)).^2).^0.5];
                    disti=[0:dist0(end)/50:dist0(end)]; 
                    xn=interp1(dist0,NOUR.xnour(i_n,:),disti);
                    yn=interp1(dist0,NOUR.ynour(i_n,:),disti);
    
                    idxm=[];
                    distmin=[];
                    dxnour0=[];
                    dsmc={};
                    for i_mc=1:COAST.n_mc
                        [ xc,yc,~,~,~ ] = get_one_polygon( COAST.x_mc,COAST.y_mc,i_mc );
                        idxn=[];
                        distn=[];
                        for kk=1:length(xn)
                        distm=((xc-xn(i_n,kk)).^2+(yc-yn(i_n,kk)).^2).^0.5;
                        distn(kk)=min(distm);
                        idxn(kk)=find(distm==min(distm),1);
                        end
                        idxm(i_mc,:)=idxn;
                        distmin(i_mc,:)=[distn]; 

                        % Compute grid cell size
                        s=[0,cumsum((diff(xc).^2+diff(yc).^2).^0.5)];
                        ds0=diff(s);
                        ds=[ds0(1),(ds0(1:end-1)+ds0(2:end))/2,ds0(end)]; 
                        % Check if coastal section is cyclic
                        cyclic=get_cyclic(xc,yc,COAST.ds0);
                        if cyclic
                            ds(1)=(ds0(1)+ds0(end))/2;
                            ds(end)=(ds0(1)+ds0(end))/2; 
                        end
                        dsmc{i_mc}=ds;
                    end
                    [dmin,idmin]=min(distmin,[],1);
                    ids0=find(diff(idmin)~=0);
                    ids1=[1,sort([ids0,ids0+1]),length(idmin)]; % find the start and end points
                    idn20=reshape(ids1,[2,length(ids1)/2])';
                    idele=idmin(ids1(1:2:end))';
                    idn=cell(COAST.n_mc,1);
                    idx=cell(COAST.n_mc,1);
                    dxnour=zeros(COAST.n_mc,1);
                    for nn=1:length(idele)
                        idn{idele(nn)}=[idn{idele(nn)},[idn20(nn,1):idn20(nn,2)]];
                        idx{idele(nn)}=[idx{idele(nn)},[min(idxm(idele(nn),[idn20(nn,1):idn20(nn,2)])):max(idxm(idele(nn),[idn20(nn,1):idn20(nn,2)]))]];
                    end
                    % determine the total length of coastline that is nourished per coastal element
                    for i_mc=1:length(idx)
                        dxnour(i_mc)=dxnour(i_mc)+sum(dsmc{i_mc}(idx{i_mc}));
                    end
                end

                % add nourishment 
                % nour [m3/m/year]
                % make sure to use only the fraction that is nourished at the considered coastal segment i_mc
                fraction = dxnour./max(sum(dxnour),1e-5);
                if ~isempty(idx{COAST.i_mc}) || sum(dxnour)==0
                    NOUR.rate_density(idx{COAST.i_mc})=NOUR.rate_density(idx{COAST.i_mc})+NOUR.rate(i_n)./dxnour(COAST.i_mc).*fraction(COAST.i_mc);
                end
            end
        end
    end
end  

function [ COAST] = merge_coastlines_mc(COAST,varargin)
% function [ COAST ] = merge_coastlines_mc(COAST,varargin)
%
% merges multiple segments of coastline (islands)
%
% INPUT:
%    COAST
%         .x_mc   : x-coordinates of the coastline elements 
%         .y_mc   : y-coordinates of the coastline elements 
%         .ds0    : grid cell size [m]
% 
% OUTPUT:
%    COAST
%         .x_mc   : x-coordinates of the coastline elements (after splitting/merging)
%         .y_mc   : y-coordinates of the coastline elements (after splitting/merging)
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
%   License along with this library. If not, see <http://www.gnu.org/licenses>
%   --------------------------------------------------------------------

    % debug info
    x_mc_orig=COAST.x_mc;
    y_mc_orig=COAST.y_mc;
    COAST.mergegrid=0;

    %% Debug information
    debug=0; % writes a memory dump to file for exceptional cases
    debugplot=0;
    sameelement=0;
    if nargin==2
        if ischar(varargin{1})
            sameelement=1;
        else
            debugplot=varargin{1};
        end
    end
    
    %% Determine epsilon and grid size
    eps=0.1;
    [ ~,~,n_mc,~,~ ] = get_one_polygon( COAST.x_mc,COAST.y_mc,1 );
    ds0=COAST.ds0;
    if ~isscalar(ds0)
        ds0=min(ds0(:,3));
    end

    for i_mc=1:n_mc-1+sameelement
        for j_mc=i_mc+1-sameelement:n_mc
            %fprintf('%3.0f %3.0f %3.0f \n',i_mc,j_mc,length(COAST.x_mc));
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% get elements i and j
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ xi,yi,~,~,~ ] = get_one_polygon( COAST.x_mc,COAST.y_mc,i_mc );
            [ xj,yj,~,~,~ ] = get_one_polygon( COAST.x_mc,COAST.y_mc,j_mc );
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Check if they are the same
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            same_element=0;
            if i_mc==j_mc || (abs(sum(xi)-sum(xj))<1e-6 && abs(sum(yi)-sum(yj))<1e-6)
                same_element=1;
            end
            
            if length(xi)>2 && length(xj)>2 % && isempty(IDS)

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% determine if elements are cylical or open
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                cyclici=hypot(xi(end)-xi(1),yi(end)-yi(1))<eps;
                cyclicj=hypot(xj(end)-xj(1),yj(end)-yj(1))<eps;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% determine if elements are clockwise
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                cwi=1;
                if cyclici
                    cwi=get_clockpoly(xi,yi);
                end
                cwj=1;
                if cyclicj
                    cwj=get_clockpoly(xj,yj);
                end
                
                % close elements if needed by making sure the last and first point are the same
                if cwi<1 && xi(end)~=xi(1) && yi(end)~=yi(1)
                    xi(end)=xi(1);
                    yi(end)=yi(1);
                end
                if cwj<1 && xj(end)~=xj(1) && yj(end)~=yj(1)
                    xj(end)=xj(1);
                    yj(end)=yj(1);
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% switch element i and j, if element i is cyclical while j is open
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if cyclici&&~cyclicj
                    % switch i and j, saves one case
                    temp=xi;xi=xj;xj=temp;
                    temp=yi;yi=yj;yj=temp;
                    temp=cyclici;cyclici=cyclicj;cyclicj=temp;
                    temp=cwi;cwi=cwj;cwj=temp;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% compute cumulative distances si,sj
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                si=get_cumdist(xi,yi); % <- not used?
                sj=get_cumdist(xj,yj); % <- not used?
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% find intersections between the two sections
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [xx,yy,indc,inds,indi,indj,ui,uj]=get_intersections(xi,yi,xj,yj);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Clean up crossings if two the same polygons are used
                % remove elements that are exactly the same and not internal 
                % crossings if 2 the same elements are used
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if same_element==1 
                    eps0=1e-6; 
                    IDdiff=abs(indc-inds)>eps0 | (ui>eps0 & ui<1-eps0); 
                    xx=xx(IDdiff); 
                    yy=yy(IDdiff); 
                    indc=indc(IDdiff); 
                    inds=inds(IDdiff); 
                    indi=indi(IDdiff); 
                    indj=indj(IDdiff); 
                    ui=ui(IDdiff); 
                    uj=uj(IDdiff); 
                end 

                IDS=[]; % elements used afterwards
                if length(xx)>1 && length(xx)~=3 && cwi==1
                    COAST.mergegrid=1;
                    % fprintf('merging elements %1.0f and %1.0f \n',i_mc,j_mc);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% Find all the crossings of the polygon (xcr, ycr)
                    %  as well as the indices of the vertices (i.e. from indi to indi+1)
                    %  including the fraction of length of the vertex (ui and uj, for i and j)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%				
                    P=[xx;yy];
                    ipi=[1:length(indi)];
                    
                    % figure;nn=3;plot(xi(indi(nn)+[0,1]),yi(indi(nn)+[0,1]),'k.-');hold on;plot(P(1,nn),P(2,nn),'r*');axis equal                
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% Add the inter-section points
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    xi2=xi;
                    yi2=yi;
                    xj2=xj;
                    yj2=yj;
                    [idi,idi2]=sort(indi+ui);
                    [idj,idj2]=sort(indj+uj);
                    rr=0;
                    for i=idi2
                        xi2=[xi2(1:indi(i)+rr),xx(i),xi2(indi(i)+1+rr:end)];
                        yi2=[yi2(1:indi(i)+rr),yy(i),yi2(indi(i)+1+rr:end)];
                        rr=rr+1;
                        indi(i)=indi(i)+rr;% max(rr-1,0);
                    end
                    rr=0;
                    for j=idj2
                        xj2=[xj2(1:indj(j)+rr),xx(j),xj2(indj(j)+1+rr:end)];
                        yj2=[yj2(1:indj(j)+rr),yy(j),yj2(indj(j)+1+rr:end)];
                        rr=rr+1;
                        indj(j)=indj(j)+rr;% max(rr-1,0);
                    end
                    
                    %if cyclici==1  || same_element==1 
                    [indi,IDsorted]=sort(indi);
                    indj=indj(IDsorted);
                    P=P(:,IDsorted);
                    xx=xx(IDsorted);
                    yy=yy(IDsorted);
                    %end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% GET SEPERATE SECTIONS OF COAST AND COMBINE
                    % THe following rules are obeyed:
                    % 1) Start with a 'clockwise' polygon (kk=0 or kk=1) from each intersection point (P) separately
                    % 2) Alternatively use the i and j polygon -> switch at intersection points (P)
                    % 3) Always go in positive direction along the paths of the polygons
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    xnew2={}; 
                    ynew2={};
                    pp2=0;
                    ids2=zeros(length(indi));
                    pp=length(indi);
                    for pp=1:length(indi)
                        %figure(15);clf;plot(xi2,yi2,'k.-');hold on;plot(xj2,yj2,'r.-');
                        %plot(P(1,:),P(2,:),'ks');

                        id1=indi(pp);
                        id2=find(indi(indi~=id1)>id1);

                        xi=xi2;
                        yi=yi2;
                        xj=xj2;
                        yj=yj2;
                        id1 = ipi(pp); 
                        ids = [];
                        % start with xi,yi if it is CW+, otherwise start with xj,yj
                        kk=1;
                        kkstart=nan;
                        if cyclicj==0 & cyclici==1
                            kk=0;
                        end
                        xnew=[];
                        ynew=[];
                        circle=0;
                        ff=0;
                        while ~circle
                            % add new points to k-element, alternating between input from element i and j
                            % with the alternation of the elements taking place at each/next crossing point 
                            if mod(kk,2)==1
                                indk0 = indi(id1);
                                indk = indi;
                                xk = xi;
                                yk = yi;
                                cwk=cwi;
                                cyclick=cyclici;
                            else
                                indk0 = indj(id1);
                                indk = indj;
                                xk = xj;
                                yk = yj;
                                cwk=cwj;
                                cyclick=cyclicj;
                            end
                            
                            %fprintf('using element %1.0f (indk0=%1.0f)\n',mod(kk,2)==1,indk0)
                            if isempty(ids)
                                id1start=id1;
                            end
                            
                            % find the elements id3 that need to be used of the k-element
                            if ~isempty(indk)
                                indklarger=indk(indk>indk0);
                                indksmaller=indk(indk<indk0);
                                % plot(xk(indk0),yk(indk0),'k*')
                                if ~isempty(indklarger)
                                    % use next line section of x,y (in positive direction) until next intersection point of the k-element
                                    id2 = find(indk==min(indklarger),1);
                                    id3 = [indk0:indk(id2)];
                                    if id2<=length(xx)
                                        xnew = [xnew, xk(id3), xx(id2)];
                                        ynew = [ynew, yk(id3), yy(id2)];
                                    else
                                        xnew = [xnew, xk(id3)];
                                        ynew = [ynew, yk(id3)];
                                    end
                                    %plot(xnew,ynew,'b-');
                                    %plot(xk(id3),yk(id3),'c-');
                                    %plot(P(1,id2),P(2,id2),'m*');hold on
                                    %plot(xk(indk0),yk(indk0),'ko')
                                    %plot([xk(id3), P(1,id2)],[yk(id3), P(2,id2)],'m-');hold on;
                                elseif ~isempty(indksmaller)
                                    % this is a special case for the last part of the k-element of an open or cylical element
                                    % use last line section of x,y & first section of x,y 
                                    % until first intersection point of the k-element, with lowest index
                                    id2 = find(indk==min(indksmaller),1);
                                    id3a = [indk0+1:length(xk)];
                                    id3b = [1:indk(id2)];
                                    if ~cyclick & mod(kk,2)==kkstart
                                        xnew = [xnew, xk(id3a)];
                                        ynew = [ynew, yk(id3a)];
                                    elseif ~cyclick 
                                        xnew = [xnew, xk(id3a), nan, xk(id3b), xx(id2)];
                                        ynew = [ynew, yk(id3a), nan, yk(id3b), yy(id2)];
                                    else
                                        xnew = [xnew, xk(id3a), xk(id3b), xx(id2)];
                                        ynew = [ynew, yk(id3a), yk(id3b), yy(id2)];
                                        % xnew = [xnew, xk(id3b), xx(id2), xk(id3a)];
                                        % ynew = [ynew, yk(id3b), yy(id2), yk(id3a)];
                                    end
                                    %plot([xk(id3a), xk(id3b), P(1,id2)],[yk(id3a), yk(id3b), P(2,id2)],'y-');hold on;
                                    %plot(xnew,ynew,'y-');hold on;
                                    %plot(P(1,id2),P(2,id2),'m*')
                                    %plot(xk(indk0),yk(indk0),'ko')    
                                else
                                    fprintf('Special case for merging (e.g. with only 1 intersection point)\n')
                                    indi=[];   % wait until it is not exactly on a point / or on the line
                                    id2=1;
                                end
                            end
                            
                            % use next intersection point
                            id1=mod(id1,length(indi))+1; %id2;
                            % switch the i and j polygon that is used
                            kk=kk+1;

                            if ~isempty(find(ids2(ids)==2)) || isempty(xnew)
                                xnew=-1e10;
                                ynew=-1e10;
                                circle=1;
                            end
                            % stopping criterion
                            %fprintf('%8.3f %8.3f %8.3f\n',pp,xnew(1)-xnew(end),ynew(1)-ynew(end))
                            if isempty(id2) & ~isempty(find(ids==id2)) ...                  % if a intersection point has been used already 2 times
                                      || ((xnew(1)-xnew(end)).^2+(ynew(1)-ynew(end)).^2).^0.5<eps ... % if begginning and end point are the approximately the same
                                      || ((xnew(end)-xj(end)).^2+(ynew(end)-yj(end)).^2).^0.5<eps ... % if the xnew,ynew end-point corresponds to end point of the original input line xj,yj
                                      || ((xnew(end)-xi(end)).^2+(ynew(end)-yi(end)).^2).^0.5<eps ... % if the xnew,ynew end-point corresponds to end point of the original input line xi,yi
                                      || (((xnew(end)-xx(id1start)).^2+(ynew(end)-yy(id1start)).^2).^0.5<eps && (cyclici || cyclicj)) 
                                      %|| (xnew(1)==xnew(end) && ynew(1)==ynew(end)) 
                                circle=1;
                                ids2(ids)=ids2(ids)+1;
                            else
                                ids  = [ids,id2];
                            end

                            % debug plot
                            %debugplot=1;
                            if debugplot==1
                                figure(pp);clf;plot(xi,yi);hold on;plot(xj,yj,'r-');plot(xj(1),yj(1),'go');plot(xj(2),yj(2),'mo');plot(xi(1),yi(1),'go');plot(xi(2),yi(2),'mo');
                                plot(xnew,ynew,'c--');plot(xnew,ynew,'c.');
                                plot(xx(pp),yy(pp),'ks')                            
                            end
                            ff=ff+1;

                            % make sure that closed polygons are recognized even if the initial eps was too strict  
                            % re-do the concatenation if with a larger eps if it is unsuccessful 
                            if ff>50
                                eps=2*eps;
                                cyclici=hypot(xi(end)-xi(1),yi(end)-yi(1))<eps;
                                cyclicj=hypot(xj(end)-xj(1),yj(end)-yj(1))<eps;
                                xi=xi2;
                                yi=yi2;
                                xj=xj2;
                                yj=yj2;
                                id1 = ipi(pp); 
                                ids = [];
                                % start with xi,yi if it is CW+, otherwise start with xj,yj
                                kk=1;
                                kkstart=nan;
                                if cyclicj==0 & cyclici==1
                                    kk=0;
                                end
                                xnew=[];
                                ynew=[];
                                circle=0;
                                ff=0;
                            end
                        end

                        % make sure to combine the first and second half of a not-closed polygon ('line'), such that it is one line again
                        idnan=find(isnan(xnew),1,'last');
                        if cyclici==0 && cyclicj==0 & ~isempty(idnan)
                            % no action on nans if 2 open line elements are used
                        else                               
                            % concatenate line sections, where the nan seperation indicates the beginning or end of section
                            getnans=1;
                            while getnans==1
                                idnan=find(isnan(xnew),1,'last');
                                if isempty(idnan)
                                    getnans=0;
                                else
                                    xnew=[xnew(idnan+1:end),xnew(1:idnan-1)];
                                    ynew=[ynew(idnan+1:end),ynew(1:idnan-1)];
                            end
                            end

                            % close polygon if it is cyclic
                            if cyclick==1 & ~(cyclici==0 || cyclicj==0)
                                xnew=[xnew,xnew(1)];
                                ynew=[ynew,ynew(1)];
                            end
                        end

                        % remove any sequential coordinates that are the same
                        uniquepoints=[1,(diff(xnew).^2+diff(ynew).^2)]~=0;
                        xnew=xnew(uniquepoints);
                        ynew=ynew(uniquepoints);

                        % get new output x,y sections
                        % split open elements at remaining location of nan
                        idnan=find(isnan(xnew));
                        idele=[0,idnan,length(xnew)+1];
                        for zz=1:length(idnan)+1
                            if xnew(1)~=-1e10
                            pp2=pp2+1;
                            xnew2{pp2}=xnew(idele(zz)+1:idele(zz+1)-1);
                            ynew2{pp2}=ynew(idele(zz)+1:idele(zz+1)-1);  
                            end
                        end
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% determine which sections to use
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % detect clockwise and anti-clockwise
                    IDclock=[];
                    for gg=1:length(xnew2)
                        % turning in CW
                        IDclock(gg) = get_clockpoly(xnew2{gg},ynew2{gg});
                    end
                    IDcw=find(IDclock>=0);
                    IDccw=find(IDclock<0);

                    % check if elements are cyclic or open
                    IDopen=[];
                    IDcyclic=[];
                    for pp=1:length(xnew2)
                        IDcyclic(pp)=1;
                        if length(xnew2{pp})>1
                            IDcyclic(pp)=hypot(xnew2{pp}(end)-xnew2{pp}(1),ynew2{pp}(end)-ynew2{pp}(1))<eps;
                        end
                        IDopen=find(IDcyclic==0);    
                    end

                    % differentiate for situations with or wiithout open elements                        
                    % if ~isempty(IDopen)
                    %     %% situation with two line elements
                    %     % element 1 : line element
                    %     % element 2 : line element
                    %     IDS=IDopen;
                    % 
                    %     %% situation with a line and island element
                    %     % -> only relevant when only 1 open line segment (before and after)
                    %     %    and where lake is enclosed by an island (e.g. island blocking of the entrance of a lagoon -> rare event)
                    %     % -> but not taking the lake into account for situations where a lake was already present initially, since it was opened then (i.e. when cwi or cwj was -1).
                    %     if length(IDopen)==1 && ~isempty(IDccw) && cwi~=-1 && cwj~=-1 && (cyclici==0 || cyclicj==0) && (cyclici==1 || cyclicj==1) 
                    % 
                    %         IDS=[IDopen,setdiff(IDccw,IDopen)];
                    %     elseif length(IDopen)==1 && ~isempty(IDccw) && cwi~=-1 && cyclici==0 && (cwj==-1 && cyclicj==1)
                    %         IDS=[IDopen,setdiff(IDcw,IDopen)];
                    %     elseif length(IDopen)==1 && ~isempty(IDccw) && cwj~=-1 && cyclicj==0 && (cwi==-1 && cyclici==1)
                    %         IDS=[IDopen,setdiff(IDcw,IDopen)];
                    %     elseif length(IDopen)==1 && cyclici==0 && (cwj==-1 && cyclicj==1)
                    %         % element 1 : line element
                    %         % element 2 : lake crossing coastline 2x -> creating island in-between
                    %         IDS=[IDopen,setdiff(IDcw,IDopen)];
                    %     end
                    % 
                    % else
                        %% differentiate elements in case of closed islands

                        % detect double polygons! (inpoly==1)                       
                        inpoly1=[];
                        for pp=1:length(xnew2)
                        for gg=1:length(xnew2)
                                % check if each of the polygons is inside another polygon?
                                if length(unique(xnew2{gg}+ynew2{gg}))<3
                                    inpoly1(pp,gg)=0;
                                elseif pp~=gg
                                    inpoly1(pp,gg)=sum(inpolygon(xnew2{pp},ynew2{pp},xnew2{gg},ynew2{gg}))/length(xnew2{pp});
                                    if (sum(xnew2{pp})-sum(xnew2{gg})==0 && sum(ynew2{pp})-sum(ynew2{gg})==0 && pp<gg) || length(xnew2{pp})<=1 || sum(isnan(xnew2{gg}))>0
                                        % keep the first of two indentical polygons 
                                        inpoly1(pp,gg)=0;
                                    end
                                else
                                    inpoly1(pp,gg)=0;
                                end
                            end
                        end
                        IDinsidenew=find(max(inpoly1,[],2)>0.9999999);

                        % determine if it is within one of the original polygons
                        inpoly2=[];
                        for pp=1:length(xnew2)
                            for gg=1:2
                                if gg==1
                                    xk=xi;
                                    yk=yi;
                                    cwk=cwi;
                                else
                                    xk=xj;
                                    yk=yj;
                                    cwk=cwj;
                                end
                                % check if each of the polygons is inside another polygon?
                                warning off
                                inpoly2(pp,gg)=sum(inpolygon(xnew2{pp},ynew2{pp},xk,yk))/length(xnew2{pp});
                                warning on
                            end
                        end
                        IDinsideorignal=find(max(inpoly2,[],2)>0.9999999);
                        IDinside_i=find(max(inpoly2(:,1),[],2)>0.9999999);
                        IDinside_j=find(max(inpoly2(:,2),[],2)>0.9999999);

                        % remove non-unique coastline point arrays
                        for pp=1:length(xnew2)
                            IDunique(pp,1)=1;
                            if length(unique(xnew2{pp}))*length(unique(ynew2{pp}))==1
                                IDunique(pp,1)=0;
                            end
                        end

                        % OUTERCONTOUR (not inside other contour)
                        if cwi>=0 && cwj>=0  
                            IDouter = setdiff(IDcw, IDinsideorignal);
                            if isempty(IDouter) && ~isempty(indi)
                                L1L2 = [length(xnew2{1}),length(xnew2{2})];
                                IDouter = setdiff(IDcw, find(L1L2==min(L1L2)));
                            elseif isempty(IDouter)
                                IDouter = IDcw;
                            end
                        else
                            %if cwi==-1 || cwj==-1
                            IDouter = IDcw;
                        end

                        % INNERCONTOUR turning in CCW : LAKES (inner contours) 
                        if cwi>=0 && cwj>=0
                            IDinner = IDccw;
                        elseif cwi>=0 && cwj==-1
                            IDinner = []; %setdiff(IDccw,IDinside_j);
                        elseif cwi==-1 && cwj>=0
                            IDinner = []; %intersect(IDccw,find(max(inpoly,[],2)>=0.9999999&min(inpoly,[],2)>-0.999999),IDccw);
                        elseif cwi==-1 && cwj==-1
                            IDinner = setdiff(IDccw,IDinsidenew);
                        end
                        IDS=[IDouter(:);IDinner(:)];
                        %figure;plot(COAST.x_mc,COAST.y_mc)

                        % filter out unnecessary polygons if there is a double intersection of a single small (cyclic) polygon with an unclosed (non-cyclic) polygon
                        if length(IDS)>2
                            if length(IDinner)>1
                                len1=[];
                                for ss=1:length(IDinner)
                                   len1(ss)=length(xnew2{IDinner(ss)});
                                end
                                IDinner=IDinner(find(len1==max(len1),1)); 
                            end
                            if length(IDouter)>=2
                                len1=[];
                                for ss=1:length(IDouter)
                                   len1(ss)=hypot(xnew2{IDouter(ss)}(end)-xnew2{IDouter(ss)}(1),ynew2{IDouter(ss)}(end)-ynew2{IDouter(ss)}(1));
                                end
                                IDouter=IDouter(find(len1==max(len1),1));
                            end
                            IDS=[IDouter(:);IDinner(:)];
                        end
                    % end
                    IDS=unique(IDS);
                    if ~isempty(IDopen)
                        IDS=[IDopen;setdiff(IDS,IDopen)];
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% DEBUG PLOTS
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    debugplot=0;
                    if debugplot==1
                        figure(30);clf
                        plot(xi,yi,'.-b',xj,yj,'.-r',xx,yy,'ok');
                        hold on

                        for gg=1:length(IDouter)
                            plot(xnew2{IDouter(gg)},ynew2{IDouter(gg)},'y-','LineWidth',1);
                        end
                        for gg=1:length(IDinner)
                            plot(xnew2{IDinner(gg)},ynew2{IDinner(gg)},'c-','LineWidth',1);
                        end
                        plotnumbers=0;
                        if plotnumbers==1
                            num=[1:length(xi)];
                            for i=1:length(xi)
                                text(xi(i),yi(i),num2str(num(i)));
                            end
                            num=[1:length(xj)];
                            for i=1:length(xj)
                                text(xj(i),yj(i),num2str(num(i)));
                            end
                        end
                        hold off
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% CREATE OUTPUT VARIABLES
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    xnewi=-1e10;
                    ynewi=-1e10;
                    xnewj=-1e10;
                    ynewj=-1e10;
                    if isempty(indi)
                        % nothing happened
                        xnewi=xi;
                        ynewi=yi;
                        xnewj=xj;
                        ynewj=yj;
                    elseif length(IDS)==1
                        xnewi=xnew2{IDS(1)};
                        ynewi=ynew2{IDS(1)};
                        xnewj=-1e10;
                        ynewj=-1e10;
                    elseif length(IDS)==2
                        if IDS(2)==IDS(1)
                            IDS(2)=mod(IDS(2),2)+1;
                        end
                        xnewi=xnew2{IDS(1)};
                        ynewi=ynew2{IDS(1)};
                        xnewj=xnew2{IDS(2)};
                        ynewj=ynew2{IDS(2)};
                    else 
                        if length(xx)~=1
                            disp('exceptional case')
                            if debug==1
                            save('debug_merge_mc.mat');
                            end
                        end
                        % do nothing
                        xnewi=xi;
                        ynewi=yi;
                        xnewj=xj;
                        ynewj=yj;
                    end
                    [COAST.x_mc,COAST.y_mc]=insert_section(xnewi,ynewi,COAST.x_mc,COAST.y_mc,i_mc);
                    [COAST.x_mc,COAST.y_mc]=insert_section(xnewj,ynewj,COAST.x_mc,COAST.y_mc,j_mc);
                end
            end
            
            % check if xi is very small (only 1 or less x values)
            if length(xi)<=2 
                xnewi=-1e10;
                ynewi=-1e10;
                [COAST.x_mc,COAST.y_mc]=insert_section(xnewi,ynewi,COAST.x_mc,COAST.y_mc,i_mc);
            end
            % check if xj is very small (only 1 or less x values)
            if length(xj)<=2 
                xnewj=-1e10;
                ynewj=-1e10;
                [COAST.x_mc,COAST.y_mc]=insert_section(xnewj,ynewj,COAST.x_mc,COAST.y_mc,j_mc);
            end
            % check if elements are the same!
            if (abs(sum(xi)-sum(xj))==0 && abs(sum(yi)-sum(yj))==0) && same_element==0  
                xnewj=-1e10;
                ynewj=-1e10;   
                [COAST.x_mc,COAST.y_mc]=insert_section(xnewj,ynewj,COAST.x_mc,COAST.y_mc,j_mc);
            end
        end
    end
    
    % remove unnecessary NaNs
    [COAST.x_mc,COAST.y_mc]=get_nansremoved(COAST.x_mc,COAST.y_mc);
    
    % check size of elements  -> sometimes incredibly small circles (bubbles) remain, which should be removed
    for i_mc=1:n_mc
        [ xi,yi,~,~,~ ] = get_one_polygon( COAST.x_mc,COAST.y_mc,i_mc );
        s=zeros(size(xi));
        s(1)=0;
        for i=2:length(xi)
            s(i)=s(i-1)+hypot(xi(i)-xi(i-1),yi(i)-yi(i-1));
        end
        if s(end)<ds0
            [COAST.x_mc,COAST.y_mc]=insert_section(-1e10,-1e10,COAST.x_mc,COAST.y_mc,i_mc);
        end
    end
    
    % remove the dummy elements that have merged into another element
    iddummy = find(COAST.x_mc==-1e10);
    for ii=length(iddummy):-1:1
        if iddummy(ii)>=length(COAST.x_mc)-1
            COAST.x_mc=[COAST.x_mc(1:iddummy(ii)-2)];
            COAST.y_mc=[COAST.y_mc(1:iddummy(ii)-2)];
        elseif iddummy(ii)==1
            COAST.x_mc=[COAST.x_mc(iddummy(ii)+2:end)];
            COAST.y_mc=[COAST.y_mc(iddummy(ii)+2:end)];
        else
            COAST.x_mc=[COAST.x_mc(1:iddummy(ii)-1),COAST.x_mc(iddummy(ii)+2:end)];
            COAST.y_mc=[COAST.y_mc(1:iddummy(ii)-1),COAST.y_mc(iddummy(ii)+2:end)];
        end
    end
end

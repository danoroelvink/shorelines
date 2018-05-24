function [QS]=Boundary_condition(QS,S,cyclic,tnow,it)
% this function applied for non cyclic shorelines 
% In case of more than one (non-cyclic) section with different...
% ... boundary the code should be adjusted 

%'FIXD2'...Neumann boundary
%'CTAN'...Constant orientation 
%'FUNC'...Dirichlet boundary 
%'FIXD'...Dirichlet (wall boundary)

%'PRDC2'...Periodic boundary Q+position   % testing phase
%'PRDC3'...Periodic boundary Q only       % testing phase

%% start point (left)
if it==0 &&strcmpi(S.boundary_condition_start,'CTAN')
    S.QS_start=QS(1);
end
if strcmpi(S.boundary_condition_start,'FIXD') && ~cyclic
    QS(1)=0; 
elseif strcmpi(S.boundary_condition_start,'FIXD2') && ~cyclic
    QS(1)=QS(2);
elseif strcmpi(S.boundary_condition_start,'FUNC') && ~cyclic
    BCraw=load(S.BCfile);warning off;
    BC=struct;
    BC.timenum=datenum([num2str(BCraw(:,1))],'yyyymmdd');
    BC.QS_strt = interpNANs(BCraw(:,2));
    BC.QSSi=interp1(BC.timenum,BC.QS_strt,tnow);
    QS(1)=interpNANs(BC.QSSi);
elseif strcmpi(S.boundary_condition_start,'CTAN') && ~cyclic
    QS(1)=S.QS_start;
elseif strcmpi(S.boundary_condition_start,'PRDC2') && ~cyclic
    QS(1)=QS(end-1);
elseif strcmpi(S.boundary_condition_start,'PRDC3') && ~cyclic
    QS(1)=QS(end-1);
           end
%% end point (right)

if it==0 &&strcmpi(S.boundary_condition_end,'CTAN')
    
    S.QS_end=QS(end);
end
if strcmpi(S.boundary_condition_end,'FIXD') && ~cyclic
    QS(end)=0;
elseif strcmpi(S.boundary_condition_end,'FIXD2') && ~cyclic
    QS(end)=QS(end-1);
elseif strcmpi(S.boundary_condition_end,'FUNC') &&  strcmpi(S.boundary_condition_start,'FUNC') && ~cyclic
    BC.QS_end=interpNANs(BCraw(:,3));
    BC.QSEi=interp1(BC.timenum,BC.QS_end,tnow);
    QS(end)=interpNANs(BC.QSEi);
elseif strcmpi(S.boundary_condition_end,'FUNC') &&  ~strcmpi(S.boundary_condition_start,'FUNC') && ~cyclic
    BCraw=load(S.BCfile);warning off;
    BC=struct;
    BC.timenum=datenum([num2str(BCraw(:,1))],'yyyymmdd');
    BC.QS_end=interpNANs(BCraw(:,3));
    BC.QSEi=interp1(BC.timenum,BC.QS_end,tnow);
    QS(end)=interpNANs(BC.QSEi);
elseif strcmpi(S.boundary_condition_end,'CTAN') && ~cyclic
    QS(end)=S.QS_end;
elseif strcmpi(S.boundary_condition_end,'PRDC2') && ~cyclic
               QS(end)=QS(2);
elseif strcmpi(S.boundary_condition_end,'PRDC3') && ~cyclic
               QS(end)=QS(2);
           
end

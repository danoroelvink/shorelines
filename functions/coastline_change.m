function [x,y,S]=coastline_change(S,s,QS,x,y,cyclic,n,tplot,tnow,tend,adt,nour,it,WVC)
dSds=zeros(size(s));
ndev=zeros(size(s));
dndt=zeros(size(s));
dx=zeros(size(s));
dy=zeros(size(s));
for i=1:n
    if cyclic
        im1=mod2(i-1,n);
        ip1=mod2(i+1,n);
    else
        im1=max(i-1,1);
        ip1=min(i+1,n+1);
    end
    
    dSds(i)=(QS(i)-QS(im1))*2/max(hypot(x(ip1)-x(im1),y(ip1)-y(im1)),S.ds0);
    %dSds(i)=(QS(i)-QS(im1))*2/hypot(x(ip1)-x(im1),y(ip1)-y(im1));
    %dSds(i)=(S(i)-S(im1))*2/...
    %   (hypot(x(ip1)-x(i),y(ip1)-y(i))+hypot(x(i)-x(im1),y(i)-y(im1)));
%     dQ(i)=(QS(i)-QS(im1));
%     dQ(end+1)=dQ(end);
    
    dndt(i)=-dSds(i)/S.d;
    ndev(i)=hypot(x(i)-.5*(x(im1)+x(ip1)),y(i)-.5*(y(im1)+y(ip1)));
end


% Bondary Conditions
if strcmpi(S.boundary_condition_start,'PRDC2') && ~cyclic
    dndt(1)=dndt(end-1);
%             elseif strcmpi(S.boundary_condition_start,'FIXD2') && ~cyclic
%                  dndt(1)=dndt(2);
end
if strcmpi(S.boundary_condition_end,'PRDC2') && ~cyclic
    dndt(end)=dndt(2);
%             elseif strcmpi(S.boundary_condition_end,'FIXD2') && ~cyclic
%                  dndt(end)=dndt(end-1);
end



if cyclic
    nend=n;
else
    nend=n+1;
end
for i=1:nend
    if cyclic
        im1=mod2(i-1,n);
        ip1=mod2(i+1,n);
    else
        im1=max(i-1,1);
        ip1=min(i+1,n+1);
    end
  
%     adtminss=adt*365*24*60
%     adt=adt*100;
    if adt <1e-05
        adt=1e-05;
    end
    if ~isempty(S.SLplot) && tplot~=0 %%for shorelines extraction
        
        if ~isempty(S.WVCfile)
            S.dt=min(S.tc*adt,(min((WVC.timenum(it+2)-tnow)/365,min((tend-tnow)/365,(tplot-tnow)/365))));
        else
            S.dt=min(S.tc*adt,min((tend-tnow)/365,(tplot-tnow)/365));
        end
        
    else
        if ~isempty(S.WVCfile)
            S.dt=min(S.tc*adt,(min((WVC.timenum(it+2)-tnow)/365,(tend-tnow)/365)));
        else
            S.dt=min(S.tc*adt,(tend-tnow)/365);
        end
    end
    
    
%     S.dt=3/48;
    
    dn=(dndt(i)+S.growth+nour(i)*S.nourrate)*S.dt;
    
    
    
    dx(i)=-dn*(y(ip1)-y(im1))/hypot(x(ip1)-x(im1),y(ip1)-y(im1));
    dy(i)= dn*(x(ip1)-x(im1))/hypot(x(ip1)-x(im1),y(ip1)-y(im1));
    %courant=S.dt/S.ds0/S.d*dQ/dn
%                  adt2(i)=S.ds0*S.d*dn/dQ(i)*S.Courant;
end
%        adtc2=abs(min(adt2));

%         if adtc2*1.00000001 < S.dt
% %             'HEEY'
%           %  pause
%         end
for i=1:nend
    x(i)=x(i)+dx(i);
    y(i)=y(i)+dy(i);
end
if cyclic
    x(n+1)=x(1);
    y(n+1)=y(1);
end
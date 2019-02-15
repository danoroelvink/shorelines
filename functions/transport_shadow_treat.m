function [xS,yS,shadowS,shadowS_h,shadow,philoc_cor,sphibr,S,hsbr]=transport_shadow_treat(x,y,x_mc,y_mc,x_hard,y_hard,phiw,philoc,sphibr,S,hsbr)

philoc_cor=philoc;

[ xS,yS,shadowS ]        = find_shadows_mc( x,y,[x_mc],[y_mc],phiw,0 );

if ~isempty(x_hard)&&~isempty(x)
    nhard=length(find(isnan(x_hard)))+1;
    for ihard=1:nhard
        [x_h,y_h]=get_one_polygon(x_hard,y_hard,ihard);
        [ ~,~,shadowS_h ]=find_shadows_mc(x,y,x_h,y_h,phiw,1);
        S.Hso(shadowS_h)=S.Hso(shadowS_h)*S.wavetransm(ihard);
        if ~isempty(hsbr)
            hsbr(shadowS_h)=hsbr(shadowS_h)*S.wavetransm(ihard);
        end
    end
    [ ~,~,shadowS_h ]      = find_shadows_mc( x,y,[x_hard],[y_hard],phiw,1 );
    shadow=zeros(size(x));
    shadow(1:end-1)=shadowS_h;
    shadow(2:end)=min(shadow(2:end),shadowS_h);
    
else
    shadowS_h=[];
    shadow=zeros(size(x));
end
philoc_cor(shadowS)=0;
%    philoc_cor(shadowS_h)=0;

sphibr(shadowS)=0;
%    sphibr(shadowS_h)=0;

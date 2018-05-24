function [ v] = get_one_polyvar( v_mc,i_mc )
nans=find(isnan(v_mc));
n_mc=length(nans)+1;
if isempty(nans)
    i1=1;
    i2=length(v_mc);
else
    if i_mc==1
        i1=1;
        i2=nans(i_mc)-1;
    elseif i_mc==n_mc;
        i1=nans(i_mc-1)+1;
        i2=length(v_mc);
    else
        i1=nans(i_mc-1)+1;
        i2=nans(i_mc)-1;
    end
end
v=v_mc(i1:i2);

end


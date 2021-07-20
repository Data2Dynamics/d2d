% $Header: svn://172.19.32.13/trunk/AMIGO2R2016/Kernel/OPT_solvers/eSS/ssm_round_int.m 770 2013-08-06 09:41:45Z attila $
function xrounded=round_int(x,index,x_L,x_U)

index=size(x,2)-index+1;

index=index:size(x,2);

for i=1:size(x,1)
    x(i,index)=x_L(index)+floor(0.5+x(i,index)-x_L(index));
end

xrounded=x;
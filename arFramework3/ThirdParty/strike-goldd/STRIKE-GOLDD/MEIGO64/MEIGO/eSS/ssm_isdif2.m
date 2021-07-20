% $Header: svn://172.19.32.13/trunk/AMIGO2R2016/Kernel/OPT_solvers/eSS/ssm_isdif2.m 770 2013-08-06 09:41:45Z attila $
function [f,ind,ind2]=isdif2(x,group,tol,flag)
f=0;
ind=[];
ind2=[];

for i=1:size(group,1)
    num=abs(x-group(i,:));
    denom=min([abs(x); abs(group(i,:))]);
    denom(abs(denom)<eps)=1;
    diference=num./denom;
    aaa=find(diference>=tol);
    
    if flag==1
        if isempty(aaa)
            f=f+1;
            ind=[ind i];
        end
    elseif flag==2
        if length(aaa)~=length(x)
            %If this happens, not all the variables comply with the
            %tolerance
            ind=[ind i];
        else
            ind2=[ind2 i];
        end
    end
end
return

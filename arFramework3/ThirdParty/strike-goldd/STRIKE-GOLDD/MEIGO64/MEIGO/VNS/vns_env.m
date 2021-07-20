function [f]=vns_env(x,y)

global cl cu ncu ncl fobj

if ~isempty(cu)
    [f,c]=feval(fobj,x);
    
    if size(c,2)==1;
        c=c';
    end
    
    temp1=c(ncu)-cu(ncu);
    temp2=cl(ncl)-c(ncl);

    temp1(temp1<0)=0;
    temp2(temp2<0)=0;

    f=f+1e6*(sum(temp1)+sum(temp2));
else
    [f]=feval(fobj,x);
end
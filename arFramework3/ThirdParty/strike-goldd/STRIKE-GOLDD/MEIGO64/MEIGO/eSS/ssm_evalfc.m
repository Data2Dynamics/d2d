% $Header: svn://172.19.32.13/trunk/AMIGO2R2016/Kernel/OPT_solvers/eSS/ssm_evalfc.m 1080 2013-11-08 12:17:11Z davidh $
function [value,value_penalty,pena,nlc,include,x] = ssm_evalfc(x,x_L,x_U,fobj,nconst,c_L,c_U,...
    tolc,weight,int_var,bin_var,varargin)

include=1;

%Transform into a row
if size(x,2)==1
    
    x=x';
    
end


if (int_var) || (bin_var)
    
    x=ssm_round_int(x,int_var+bin_var,x_L,x_U);
    
end

%Adjust to bounds
lower=find(x<x_L);
upper=find(x>x_U);

x(lower)=x_L(lower);
x(upper)=x_U(upper);

if nconst
    
    [value,nlc]=feval(fobj,x,varargin{:});
    
    %Transform into a row
    if size(nlc,2)==1
        
        nlc=nlc';
        
    end
    
    pena=ssm_penalty_function(x,nlc,c_L,c_U,tolc);
    value_penalty=value+weight*pena;
       
else
    
    value=feval(fobj,x,varargin{:});
    nlc=0;
    pena=0;
    value_penalty=value;
    
end

if isinf(value) || isnan(value)
    
    include=0;
    
end

return

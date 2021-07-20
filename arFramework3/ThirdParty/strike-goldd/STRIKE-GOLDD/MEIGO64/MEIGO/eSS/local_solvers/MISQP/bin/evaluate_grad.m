function [ df,dg,nfunc ] = evaluate_grad( x,f,g,ncont,nfunc,m,df,dg,FUN,neq,fobj,varargin)


%EVAL_GRAD Gradient evaluation
%
    eps = 1.0d-7;

    %gradients only for continuous variables
    index=ncont;
    feps=0;
    wa1=zeros(1,m);
    %nfunc=0;
    for i=1:index
        epsa = eps*max(1.0d-5,abs(x(i)));
        epsi=1.0d0/epsa;
        x(i)=x(i)+epsa;
        [ feps,wa1 ] = ...
            feval(FUN, x,neq,fobj,varargin{:});

        nfunc=nfunc+1;

        df(i)=epsi*(feps-f);
        for j=1:m
            dg(j,i)=epsi*(wa1(j)-g(j));
        end
        x(i)=x(i)-epsa;
    end
end



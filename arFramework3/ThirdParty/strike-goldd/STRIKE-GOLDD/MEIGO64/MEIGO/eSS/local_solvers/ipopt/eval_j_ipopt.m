function [J] = eval_j_ipopt(x,returnStructureOnly)
global n_fun_eval hjfun params

neq=0;
index=length(x);
%EVAL_GRAD Gradient evaluation
%
%
eps = 1.0d-7;

[ f,g ] = feval(hjfun,x,params{:} );
nconst=length(g);
n_fun_eval=n_fun_eval+1;

if returnStructureOnly
    J = sparse(ones(nconst,index));
else
    feps=0;
    wa1=zeros(1,nconst);

    dg=zeros(nconst,index);

    for i=1:index
        epsa = eps*max(1.0d-5,abs(x(i)));
        epsi=1.0d0/epsa;
        x(i)=x(i)+epsa;
        [feps,wa1] = feval(hjfun,x,params{:});
        n_fun_eval=n_fun_eval+1;

        %df(i)=epsi*(feps-f);
        for j=1:nconst
            %Vectorize this!!!!
            dg(j,i)=epsi*(wa1(j)-g(j));
        end
        x(i)=x(i)-epsa;
    end
    J=sparse(dg);
end
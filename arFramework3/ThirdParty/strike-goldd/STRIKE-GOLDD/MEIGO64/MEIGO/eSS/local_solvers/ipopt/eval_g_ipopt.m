function df = eval_g_ipopt(x)

global n_fun_eval hjfun params

neq=0;
index=length(x);
%EVAL_GRAD Gradient evaluation
%  
%    
eps = 1.0d-7;

[f] = feval(hjfun,x,params{:} );
n_fun_eval=n_fun_eval+1;
feps=0;
df=zeros(1,index);
for i=1:index
    epsa = eps*max(1.0d-5,abs(x(i)));
    epsi=1.0d0/epsa;
    x(i)=x(i)+epsa;
    [feps] = feval(hjfun,x,params{:});
    n_fun_eval=n_fun_eval+1;
                     
    df(i)=epsi*(feps-f);
%     for j=1:nconst
%         dg(j,i)=epsi*(wa1(j)-g(j));
%     end
    x(i)=x(i)-epsa;
end
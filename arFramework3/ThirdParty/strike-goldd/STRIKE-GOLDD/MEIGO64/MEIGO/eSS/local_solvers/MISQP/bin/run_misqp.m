function [x,f,nfunc,ifail,g]=run_misqp(ncont,nint,nbin,m,xl,xu,x0,acc,FUN,neq,fobj,varargin)

n=nint+ncont+nbin;

me=neq;


x=x0;

% ALLOCATION MEMORY LOCAL SOLVER MISQP
f_ls=0;
g_ls=zeros(1,m+me+1);
x_ls=zeros(1,n+1);
xl_ls=zeros(1,n+1);
xu_ls=zeros(1,n+1);
u=zeros(1,m+n+n);
b=zeros(n+1,n+1);
dg=zeros(m+me+1,n+1);
df=zeros(1,n+1);

% working arrays
lrw=[];
rw=[];
liw=[];
iw=[];
llw=[];
lw=[];


ifail=0;
maxit=500;
maxpen=50;
maxund=50;
mode=0;
%iprint=0;
%accqp = 1.0d-12;
accqp = 1e-6;
%acc=1e-3;

resopt=1;
nonmon=2;


mme= m + me;
iprint=0;

nfunc = 0;
ngrad = 0;

% Calculate function
nfunc=nfunc+1;

[ f,g ] = feval(FUN, x,neq,fobj,varargin{:} );  % here you have to specify your function
% which calculates the function
% values


[ df,dg,nfunc ] = evaluate_grad( x,f,g,ncont,nfunc,m,df,dg,FUN,neq,fobj,varargin{:});
ngrad = ngrad+1;

f_ls=double(f);
g_ls(1:m)=double(g(1:m));
x_ls(1:n)=double(x);
xl_ls(1:n)=double(xl);
xu_ls(1:n)=double(xu);

%keyboard

while(ifail <= 0)
    %   keyboard
    [x_ls,f_ls,g_ls,df,dg,u,xl_ls,xu_ls,b,ifail,rw,lrw,iw,liw,lw,llw] = ...
        misqp(m,me,n,nint,nbin,x_ls,...
        f_ls,g_ls,df,dg,u,xl_ls,xu_ls,b,acc,accqp,...
        maxit,maxpen,maxund,resopt,nonmon,iprint,mode,ifail,rw,lrw,...
        iw,liw,lw,llw);
    %  keyboard
    if (ifail == -1)
        % Calculate function
        nfunc=nfunc+1;
        [ f,g ] =...
            feval(FUN, x_ls(1:n),neq,fobj,varargin{:} );  % here you have to specify your function
        % which calculates the function
        % values

        f_ls=f;
        g_ls(1:m)=g(1:m);
    end
    if (ifail == -2)
        % Calculate gradients
        ngrad=ngrad+1;
        [df,dg,nfunc]=evaluate_grad( x_ls(1:n),f_ls,g_ls(1:m),ncont,nfunc,m,df,dg,FUN,neq,fobj,varargin{:});
        % here you have to specify your function
        % which calculates the gradients
    end
    if (ifail>=0)
        break
    end

end

x=x_ls(1:n);
f=f_ls;
g=g_ls(1:m);
%  keyboard

if size(x,1)>1
    x=x';
end

disp(ifail);

end


function [p,jh,flag,ic,lambda,hess]=solnp(FUN,pb,ib,op,lambda,hess,...
    iprint,fobj,scale,nconst,varargin)

%The Function SOLNP solves nonlinear programs in standard form:
%
%        minimize              J(P)
%        subject to            EC(P)  =0
%                   IB(:,1)<=  IC(P)  <=IB(:,2)
%                   PB(:,1)<=    P    <=PB(:,2).
%where
%
%  J       : Cost objective scalar function
%  EC      : Equality constraint vector function
%  IC      : Inequality constraint vector function
%  P       : Decision parameter vector
%  IB, PB  : lower and upper bounds for IC and P.
%
%%% USAGE
% [P,JH,FLAG,IC,LAMBDA,HESS]=SOLNP(FUN,PB,IB,OP,L,H,IPRINT,PAR1,PAR2,...)
%  Output P       : optimal decision parameters
%         JH      : optimization history (nfeval - function values)
%         FLAG    : flag indicating whether SOLNP has converged.
%         IC      : optimal values of inequality constraints
%         LAMBDA  : optimal lagrangian multipliers
%         HESSIAN : final estimate of the hessian
%
%  Input  FUN    : user-defined input function
%                   [F]=FUN(P,IT)
%                   F(1)        : cost objective function J value
%                   F(2:NEC+NIC): constraint function values o
%                   IT           : iteration number (used?)
%                                >0 - major iteration
%                                >0 - minor iteration
%                                0 - any other call
%
%         PB = [P {PB}]
%
%            P   : any initial parameter within the parameter bound
%            PB  : optional parameter bound
%                  PB(:,1) is the lower bound for P, and
%                  PB(:,2) is the upper bound for P.
%
%         IB = [{IC} IB] (optional)
%
%            IC  : (optional) best approximation values for inequality 
%                  constraints (IC) within the bounds
%            IB  : inequality constraint bound
%                  IB(:,1) is the lower bound for IC, and
%                  IB(:,2) is the upper bound for IC.
%
%         OP=[RHO,MAJIT,MINIT,MAXFN,DELTA,TOL] (all optional with defaults)
%
%           RHO  : penalty parameter (1)
%           MAJIT: maximum number of major iterations (10)
%           MINIT: maximum number of minor iterations (10)
%           MAXFN: maximum number of funcion evaluations (500*no.variables)
%           DELTA: relative step size in forward difference (1.e-5)
%           TOL  : tolerance on feasibility and optimality (1.e-4)
%
%         IPRINT : flag for printing optimization process on screen
%
% Modified by Jose-Oscar H. Sendin
% Process Engineering Group, IIM-CSIC
% Latest update: August 23, 2006
% -------------------------------------------------------------------------

% check inputs --- JOHS
if nargin < 2, error('Syntax error: SOLNP requires at least two input arguments!'); end
if nargin < 7, iprint = 0;
    if nargin < 6, hess = [];
        if nargin < 5; lambda = [];
            if nargin < 4, op = [];
                if nargin < 3, ib = []; 
                end, end, end, end, end

om=['SOLNP--> ';'         '];

% check bounds
[np,n]=size(pb);  % np = number of parameters
lpb=[1 0];
% lpb(1) = refers to the presence of parameter bounds
% lpb(2) = refers to constraints of either type
if n==1, 
    p=pb; pb=0; lpb(1)=0;
elseif n==2, 
    p=(pb(:,1)+pb(:,2))/2;
elseif n==3, 
    p=pb(:,1); pb=pb(:,2:3);
else
    error([om(1,:) 'Parameter array ''PB'' must have three columns or less']);
end

if lpb(1)==1,
    if min(pb(:,2)-pb(:,1))<=0,
        error([om(1,:) 'The lower bounds of the parameter constraints ';...
                om(2,:) 'must be strictly less than the upper bounds.  ']);
    elseif min([p-pb(:,1);pb(:,2)-p])<=0,
        error([om(1,:) 'Initial parameter values must be within the bounds']);
    end
end
% check constraints
if isempty(ib) 
    nic=0; % number of inequality constraints
    ic=0;  
else  
    [nic, n]=size(ib);
    if n==3,
        ib0=ib(:,1);ib=ib(:,[2:3]);
        if min([ib0-ib(:,1);ib(:,2)-ib0])<=0,
            error([om(1,:) 'Initial inequalities must be within the bounds']);
        end
    elseif n==2,
        if min(ib(:,2)-ib(:,1))<=0,
            error([om(1,:) 'The lower bounds of the inequality constraints';...
                   om(2,:) 'must be strictly less than the upper bounds.  ']);
       end
       ib0=(ib(:,1)+ib(:,2))/2;
   elseif n==1,
       ic=0; nic=0;
   else
       error([om(1,:) 'Inequality constraints must have 2 or 3 columns.']);
   end
   if nic>=1,
       if lpb(1)==1,
           pb=[ib; pb];
       else
           pb=ib;
       end
       p=[ib0; p]; % p contains the original variables + slack variables for
                   % the inequalities
   end
end

clear ib ib0
if ( lpb(1)==1 | nic>0 ), lpb(2)=1; end

% optimization parameters (including max fun evals)
opd=[1 10 10 1000*np 1.0e-5 1.0e-4];  % default 

if isempty(op),
    op=opd;
else
    [m,n]=size(op);
    if m>1, error([om(1,:) 'Control variables must be a vector']); end
    if n<6, op=[op(1:n) opd(n+1:6)]; end
end

op(1)=max(op(1),0); 
idx = find(op<=0);
op(idx) = opd(idx);

rho=op(1); maxit=op(2); minit=op(3); maxfn = op(4); delta=op(5); tol=op(6);

clear op opd

flag = 0;
nfeval = 0;

% evaluate initial guess

ob=feval(FUN,p(nic+1:nic+np),varargin,fobj,nconst);
nfeval = nfeval + 1;
[m,n]=size(ob);
if n>=2,
    error([om(1,:) 'Cost function must return a column vector'])
end
if m<nic+1,
    error([om(1,:) 'The number of constraints in your COST function does';...
           om(2,:) 'not match the number specified in the call to SOLNP.']);  
end 
nc=m-1; % total number of constraints
nec=nc-nic; % number of equality constraints
clear m n

j=ob(1);  % objective function value
jh=[nfeval, j]; % optimization history
t=zeros(3,1); % tolerances

% lagrange multipliers
if nc>=1,
    if ( isempty(lambda) | lambda==0 ) 
        lambda=zeros(nc,1);
    end
    constraint=ob(2:end);
    if nic>=1,
        if min([constraint(nec+1:nc)-pb(1:nic,1);...
                    pb(1:nic,2)-constraint(nec+1:nc)]) > 0,
            p(1:nic)=constraint(nec+1:nc); 
        end
        constraint(nec+1:nc)=constraint(nec+1:nc)-p(1:nic);
    end
    t(2)=norm(constraint);
    if max([t(2)-10*tol, nic])<=0, %???
        rho=0;
    end
else
    lambda = 0;
end

% hessian
if isempty(hess), 
    hess=eye(np+nic); 
end

mu=np; 
iteration=0;
% begin iteration loop
while iteration<maxit,
    iteration=iteration+1;
    op=[rho minit delta tol nec nic np lpb];
    [p,lambda,hess,mu,subnev,info]=subnp(FUN,p,op,lambda,ob,pb,hess,mu,fobj,nconst,scale,varargin{:});
    ob=feval(FUN,p(nic+1:nic+np),varargin,fobj,nconst);
    nfeval = nfeval + subnev + 1;
    
    t(1)=(j-ob(1))/max(abs(ob(1)),1);
    j=ob(1);
    if nc>=1,
        constraint=ob(2:nc+1); 
        if nic>=1,
            if min([constraint(nec+1:nc)-pb(1:nic,1);...
                        pb(1:nic,2)-constraint(nec+1:nc)]) > 0,
                p(1:nic)=constraint(nec+1:nc); 
            end
            constraint(nec+1:nc)=constraint(nec+1:nc)-p(1:nic);
        end
        t(3)=norm(constraint);
        % update penalty weight
        if t(3)<10*tol,
            rho=0; mu=min(mu, tol);
        end
        if t(3)<5*t(2), 
            rho=rho/5;
        elseif t(3)>10*t(2), 
            rho=5*max(rho, sqrt(tol)); 
        end
        if max([tol+t(1), t(2)-t(3)]) <= 0, 
            lambda=0*lambda; hess=eye(size(hess)); 
        end
        t(2)=t(3);
    end
    % optimization history
    jh=[jh; nfeval, j];
    if iprint > 0,
        fprintf(iprint,'iteration: %i; nfeval: %i; f: %g \n', iteration,nfeval,ob(1));
        if info(1)==1, 
            fprintf(iprint, '-- SOLNP -> %s\n\t%s\n\t%s\n\n', ...
                'Redundant constraints were found. Poor',...
                'intermediate results may result.  Suggest that you',...
                'remove redundant constraints and re-OPTIMIZE.'); 
        end
        if info(2)==1, 
            fprintf(iprint, '-- SOLNP -> %s\n\n', ...
                'The linearized problem has no feasible solution.'); 
        end
        if info(3)==1, 
            fprintf(iprint, '-- SOLNP -> %s\n\t%s\n\n', ...
                'Minor optimization routine did not converge in the',...
                'specified number of minor iterations.'); 
        end
    end
    % termination
    if norm([t(1) t(2)])<=tol, 
        maxit = iteration;
        flag = 1;
    end
    if nfeval > maxfn, break, end
    if info(2)==1 & abs(t(1))<=eps & t(2)>tol, flag = -1; break, end % no progress
end

if nic>=1,
    ic=p(1:nic);
end
p=p(nic+1:nic+np);

if iprint
    if flag == 1
        fprintf(iprint,'-- Completed in %i iteration \n',iteration);
    elseif flag == -1;
        fprintf(iprint,'-- No feasible solution found \n');
    elseif nfeval > maxfn,
        fprintf(iprint,'-- Exiting after maximum number of function evaluations \n');
    else
        fprintf(iprint,'-- Exiting after maximum number of iterations \n');
    end
end
return

%  VARIABLE GLOSSARY:
%-------------------------------------------------------------------------------
%OB(1)      value of the cost objective function
%CONSTRAINT:vector of constraint values
%IB:        on input, contains the inequality constraint bounds + optionally
%           the values of the inequality constraints.  Gets converted to a 
%           NIC x 2 matrix of inequality constraint bounds.
%IC:        NIC x 1 vector of inequality constraints
%ITERATION: index for major iterations
%J:         previous value of the cost objective function
%JH:        history of the cost function
%LPB (2):   vector flag which indicates the presence of parameter bounds
%           and or inequality constraints.
%             LPB(1) refers to parameter bounds, it is 0 if there are 
%                    none, 1 if there are one or more.  
%             LPB(2) refers to constraints of either type.
%M:         number of rows of a matrix, usually temporary
%N:         number of columns of a matrix, usually temporary
%NC:        total number of constraints (=NEC+NIC)
%NEC:       number of equality constraints
%NIC:       number of inequality constraints
%NP:        number of parameters
%OPD:       vector of default optimization control variables
%OP:        vector of control variables.
%             It is passed in as:
%               [RHO MAJIT MINIT DELTA TOL] (all optional)
%             It is passed to ISISUBOPT as:
%               [RHO MAJIT MINIT DELTA TOL NEC NIC NP LPB(1) LPB(2)]
%P:         On input, contains the parameters to be optimized.  During the 
%           optimization this vector contains: [PIC;P]  where PIC contains 
%           pseudo parameters corresponding to the inequality constraints.
%PB:        on input, optionally contains the parameter bounds + optionally
%           the values of the parameters (one or the other or both can be 
%           specified).  Gets converted to a NPB x 2 matrix of parameter bounds.
%T (3):     vector of computed tolerances during optimization.
%             T(1) is the difference of the objective values between two 
%                  consecutive iterations
%             T(2) is NORM(CONSTRAINT) before a major iteration
%             T(3) is NORM(CONSTRAINT) after a major iteration


%-------------------------------------------------------------------------------
%  SUBNP
%-------------------------------------------------------------------------------
function [p,y,hess,mu,nfeval,info]=subnp(FUN,p0,op,yy,ob,pb,hess,mu,fobj,nconst,scale_temp,varargin)

rho=op(1); maxit=op(2); delta=op(3); tol=op(4);   nec=op(5); 
nic=op(6); np=op(7);    nc=nec+nic;  npic=np+nic; lpb=op(8:9); ch=1;
clear op
alp=[0 0 0];
nfeval = 0;

info = zeros(1,3);

% make the scale for the cost, the equality constraints, the inequality
% constraints, and the parameters
if nec>0,
    scale=[ob(1);ones(nec,1)*max(abs(ob(2:nec+1)))];
else        
    scale=1;
end
if lpb(2)==0,
    scale=[scale; p0];
else        
    scale=[scale; ones(size(p0))];
end

%JAE 10/07/2006
% scale_temp=varargin{2};
scale=[scale_temp.f scale_temp.g scale_temp.x]';
%


scale=min(max(abs(scale),tol),1/tol);




% scale the cost, the equality constraints, the inequality constraints, 
% the parameters (inequality parameters AND actual parameters), 
% and the parameter bounds if there are any
% Also make sure the parameters are no larger than (1-tol) times their bounds
ob=ob./scale(1:nc+1);  
p0=p0./scale(nec+2:nc+np+1);
if lpb(2)==1,
    if lpb(1)==0,
        mm=nic;
    else 
        mm=npic;
    end
    pb=pb./[scale(nec+2:nec+mm+1) scale(nec+2:nec+mm+1)];
end

% scale the lagrange multipliers and the Hessian
if nc>=1,
    yy=scale(2:nc+1).*yy/scale(1);
end
hess=hess.*(scale(nec+2:nc+np+1)*scale(nec+2:nc+np+1)')/scale(1);

j=ob(1);
a=[zeros(nec,nic);-eye(nic)];
g=zeros(npic,1);

% evaluate gradients using forward differences
if nc>=1,
    constraint=ob(2:nc+1);
    for i=1:np, 
        p0(nic+i)=p0(nic+i)+delta;
        ob=feval(FUN,p0(nic+1:npic).*scale(nc+2:nc+np+1),varargin,fobj,nconst);
        ob = ob./scale(1:nc+1);
        g(nic+i)=(ob(1)-j)/delta;
        a(:,nic+i)=(ob(2:nc+1)-constraint)/delta;
        p0(nic+i)=p0(nic+i)-delta;
    end
    nfeval = nfeval + np;
    if nic>=1,
        constraint(nec+1:nec+nic)=constraint(nec+1:nec+nic)-p0(1:nic);
    end
    if cond(a)>1/eps,
        info(1) = 1; % 'Redundant constraints';
    end
    b=a*p0-constraint;
end

% see if the current vector is feasible for the linearized problem
if nc>=1,
    ch=-1;
    alp(1)=tol-max(abs(constraint));
    if alp(1)<=0,
        ch=1;
        if lpb(2)==0,
            p0=p0-a'*((a*a')\constraint);
            alp(1)=1;
        end
    end
    if alp(1)<=0,
        p0(npic+1)=1;
        a=[a, -constraint];
        c=[zeros(1,npic), 1];
        dx=ones(npic+1,1);
        go=1; 
        minit=0;
        while go>=tol 
            minit=minit+1;
            gap=[p0(1:mm,1)-pb(:,1),pb(:,2)-p0(1:mm,1)];
            gap=sort(gap')';
            dx(1:mm)=gap(:,1);
            dx(npic+1)=p0(npic+1);
            if lpb(1)==0,
                dx(mm+1:npic)=max([dx(1:mm);100])*ones(npic-mm,1);
            end
            y=(a*diag(dx))'\(dx.*c');
            v=dx.*(dx.*(c'-a'*y));
            if v(npic+1)>0,
                z=p0(npic+1)/v(npic+1);
                for i=1:mm,
                    if v(i)<0,
                        z=min(z,-(pb(i,2)-p0(i))/v(i));
                    elseif v(i)>0, 
                        z=min(z,(p0(i)-pb(i,1))/v(i)); 
                    end
                end
                if z>= p0(npic+1)/v(npic+1),
                    p0=p0-z*v; 
                else
                    p0=p0-0.9*z*v; 
                end
                go=p0(npic+1);
                if minit >= 10, 
                    go=0; 
                end
            else
                go=0;
                minit=10;
            end
        end
        if minit>=10,
            info(2) = 1; % 'Infeasible';
        end
        a=a(:,1:npic); 
        b=a*p0(1:npic);
    end
end

clear constraint c z v gap;

p=p0(1:npic); 
y=0; 
if ch>0,
    ob=feval(FUN,p(nic+1:npic).*scale(nc+2:nc+np+1),varargin,fobj,nconst);
    ob = ob./scale(1:nc+1);
    nfeval = nfeval + 1;
end
j=ob(1);

if nic>=1,
    ob(nec+2:nc+1)=ob(nec+2:nc+1)-p(1:nic);
end

if nc>=1,
    ob(2:nc+1)=ob(2:nc+1)-a*p+b;
    j=ob(1)-yy'*ob(2:nc+1)+rho*norm(ob(2:nc+1))^2;
end

minit=0; STOP = 0;
while ~STOP 
    minit=minit+1;
    if ch>0,
        for i=1:np,
            p(nic+i)=p(nic+i)+delta;
            obm=feval(FUN,p(nic+1:npic).*scale(nc+2:nc+np+1),varargin,fobj,nconst);
            obm = obm./scale(1:nc+1);
            if nic>0,
                obm(nec+2:nc+1)=obm(nec+2:nc+1)-p(1:nic);
            end
            if nc>0,
                obm(2:nc+1)=obm(2:nc+1)-a*p+b;
                obm=obm(1)-yy'*obm(2:nc+1)+rho*norm(obm(2:nc+1))^2;
            end
            g(nic+i)=(obm-j)/delta;
            p(nic+i)=p(nic+i)-delta;
        end
        nfeval = nfeval + np;
        if nic>=1,
            g(1:nic)=0*yy(nec+1:nc);
        end
    end
    if minit>1,
        yg=g-yg;
        sx=p-sx;
        sc(1)=sx'*hess*sx;
        sc(2)=sx'*yg;
        if sc(1)*sc(2)>0,
            sx=hess*sx;
            hess=hess-(sx*sx')/sc(1)+(yg*yg')/sc(2);
        end
    end
    dx=0.01*ones(npic,1);
    if lpb(2)==1,
        gap=[p(1:mm,1)-pb(:,1),pb(:,2)-p(1:mm,1)];
        gap=sort(gap')';
        gap=gap(:,1)+sqrt(eps)*ones(mm,1);
        dx(1:mm,1)=ones(mm,1)./gap;
        if lpb(1)<=0,
            dx(mm+1:npic)=min([dx(1:mm);0.01])*ones(npic-mm,1);
        end
    end
    go=-1;
    mu=mu/10;
    while go<=0,
        % using Cholesky factorization
        [c,pdummy]=chol(hess+mu*diag(dx.*dx));
        if pdummy > 0, p0 = p; break, end % ????????????????????
        c=inv(c);
        yg=c'*g;
        if nc<=0,
            u=-c*yg;
        else 
            y=(c'*a')\yg;
            u=-c*(yg-(c'*a')*y);
        end
                
        p0=u(1:npic)+p;
        if lpb(2)==0,
            go=1;
        else
            go=min([p0(1:mm)-pb(:,1);pb(:,2)-p0(1:mm)]);
            mu=3*mu;
        end
    end
    alp(1)=0;ob1=ob;ob2=ob1;sob(1)=j;sob(2)=j;
    pt(:,1:2)=[p p];alp(3)=1.0;pt(:,3)=p0;
    ob3=feval(FUN,pt(nic+1:npic,3).*scale(nc+2:nc+np+1),varargin,fobj,nconst);
    ob3 = ob3./scale(1:nc+1);
    nfeval = nfeval + 1;
    sob(3)=ob3(1);
    if nic>=1,
        ob3(nec+2:nc+1)=ob3(nec+2:nc+1)-pt(1:nic,3);
    end
    if nc>=1,
        ob3(2:nc+1)=ob3(2:nc+1)-a*pt(:,3)+b;
        sob(3)=ob3(1)-yy'*ob3(2:nc+1)+rho*norm(ob3(2:nc+1))^2;
    end
    go=1;
    while go>tol,
        alp(2)=(alp(1)+alp(3))/2;
        pt(:,2)=(1-alp(2))*p+alp(2)*p0;
        ob2=feval(FUN,pt(nic+1:npic,2).*scale(nc+2:nc+np+1),varargin,fobj,nconst);
        ob2 = ob2./scale(1:nc+1);
        nfeval = nfeval + 1;
        sob(2)=ob2(1);
        if nic>=1,
            ob2(nec+2:nc+1)=ob2(nec+2:nc+1)-pt(1:nic,2);
        end
        if nc>=1,
            ob2(2:nc+1)=ob2(2:nc+1)-a*pt(:,2)+b;
            sob(2)=ob2(1)-yy'*ob2(2:nc+1)+rho*norm(ob2(2:nc+1))^2;
        end
        obm=max(sob);
        if obm<j,
            obn=min(sob);
            go=tol*(obm-obn)/(j-obm);
        end
        if sob(2)>=sob(1),
            sob(3)=sob(2);ob3=ob2;alp(3)=alp(2);pt(:,3)=pt(:,2);
        elseif sob(1)<=sob(3),
            sob(3)=sob(2);ob3=ob2;alp(3)=alp(2);pt(:,3)=pt(:,2);
        else
            sob(1)=sob(2);ob1=ob2;alp(1)=alp(2);pt(:,1)=pt(:,2);
        end
        if go>=tol,
            go=alp(3)-alp(1);
        end
    end
    sx=p;yg=g;ch=1;
    obn=min(sob);
    
    if j<=obn, STOP = 1; end

    reduce=(j-obn)/(1+abs(j));
    if reduce<tol, STOP = 1; end
    
    if sob(1)<sob(2),
        j=sob(1);p=pt(:,1);ob=ob1;
    elseif sob(3)<sob(2),
        j=sob(3);p=pt(:,3);ob=ob3;
    else 
        j=sob(2);p=pt(:,2);ob=ob2;
    end
    if minit >= maxit, STOP = 1; end
    clear ob1 ob2 ob3 pt;
end

p=p.*scale(nec+2:nc+np+1);  % unscale the parameter vector
if nc>=1,
    y=scale(1)*y./scale(2:nc+1); % unscale the lagrange multipliers
end
hess=scale(1)*hess./(scale(nec+2:nc+np+1)*scale(nec+2:nc+np+1)');

if reduce>tol,
    info(3) = 1; % 'Max number of minor iterations';
end
return

%  VARIABLE GLOSSARY:
%------------------------------------------------------------------------------
%
%ALP (4):   vector of ??. Alp(4) is "alpha"
%CH:        "change"
%CONSTRAINT:vector of constraint values
%G:         gradient
%GAP (2):   lower and upper "gap"
%LPB (2):   vector flag which indicates the presence of parameter bounds
%           and or inequality constraints.
%             LPB(1) refers to parameter bounds, it is 0 if there are 
%                    none, 1 if there are one or more.  
%             LPB(2) refers to constraints of either type.
%NC:        total number of constraints (=NEC+NIC)
%NEC:       number of equality constraints
%NIC:       number of inequality constraints
%NP:        number of parameters
%P0:        parameter vector on input
%P:         parameter vector on output
%SOB (3):   vector of ??
%Z (5):     vector of ??
 function [x,f,g,df,dg,u,xl,xu,b,ifail,rw,lrw,iw,liw,lw,...
 llw] = misqp(m,me,n,nint,nbin,x,...
 f,g,df,dg,u,xl,xu,b,acc,accqp,maxit,maxpen,maxund,resopt,nonmon,iprint,mode,...
 ifail,rw,lrw,iw,liw,lw,llw)
% ***********************************************************************
% 
%
%    In order to speed up the interface, this version is limited to
%    n <= 1000
%    m <= 1000
%    
%    
% 
%         An Implementation of a Trust Region Method for Solving
%             Mixed-Integer Nonlinear Optimization Problems 
% 
% 
%    MISQP solves the mixed-integer nonlinear program (MINLP)
% 
%              minimize        F(X,Y)
%              subject to      G(J)(X,Y)   =  0  , J=1,...,ME
%                              G(J)(X,Y)  >=  0  , J=ME+1,...,M
%                              XL  <=  X  <=  XU
%                              YL  <=  Y  <=  YU
% 
%    where X is a real and y an integer variable vector.
% 
%    The Fortran subroutine is an implementation of a modified sequential 
%    quadratic programming (SQP) method. Under the assumption that integer 
%    variables have a 'smooth' influence on the model functions, i.e., that 
%    function values do not change too drastically when in- or decrementing 
%    an integer value, successive quadratic approximations are applied.
%    The algorithm is stabilized by a trust region method with Yuan's second 
%    order corrections. 
% 
%    It is not assumed that the mixed-integer program is relaxable. In other 
%    words, function values are required only at integer points. The Hessian 
%    of the Lagrangian function is approximated by BFGS updates subject to 
%    the continuous variables. Derivative information subject to the integer 
%    variables is obtained by a difference formula evaluated at grid points.
% 
% 
%    USAGE:
% 
%       [x,f,g,df,dg,u,xl,xu,b,ifail,rw,iw,lw] = misqp(m,me,n,nint,nbin,x,...
%                      f,g,df,dg,u,xl,xu,b,acc,accqp,maxit,maxpen,maxund,...
%                      resopt,nonmon,iprint,mode,ifail,rw,lrw,iw,liw,lw,llw)
% 
% 
%    ARGUMENTS:
% 
%    M :       Total number of constraints.
%    ME :      Number of equality constraints.
%    N :       Number of optimization variables, continuous and integer ones.
%    NINT :    Number of integer variables, must be less than or equal to N.
%    NBIN :    Number of binary variables, must be less than or equal to N.
%    MNN :     Must be equal to M+N+N, dimensioning parameter of multiplier 
%              vector.
%    X(N+1) :    Initially, X has to contain starting values. On return, X is
%              replaced by the current iterate. In the driving program
%              the row dimension of X has to be equal to NMAX.
%              X contains first the continuous, then the integer followed by
%              the binary variables.
%    F :       F contains the actual objective function value evaluated at X.
%    G(M+ME+1) :  G contains the actual constraint function values at X.
%    DF(N+1) :  DF contains the gradient values of the objective function.
%    DG(M+ME+1,N+1) :  DG contains the gradient values of constraints
%              in the first M rows. In the driving program, the row dimension
%              of DG has to be equal to MMAX.
%    U(M+N+N) :  On return, the first M locations contain
%              the multipliers of the M nonlinear constraints, the subsequent
%              N locations the multipliers subject to the lower bounds, and the
%              final N locations the multipliers subject to the upper bounds.
%              At an optimal solution, all multipliers with respect to
%              inequality constraints should be nonnegative.
%    XL(N+1),XU(N+1) :  On input, the one-dimensional arrays XL and XU must
%              contain the upper and lower bounds of the variables, first
%              for the continuous, then for the integer and subsequently for 
%              the binary variables.
%    B(N+1,N+1) :  On return, B contains the last computed approximation
%              of the Hessian matrix of the Lagrangian function.
%              In the driving program, the row dimension of C has to be equal
%              to NMAX. 
%    ACC :     The user has to specify the desired final accuracy
%              (e.g. 1.0D-7). The termination accuracy should not be smaller
%              than the accuracy by which gradients are computed. If ACC is
%              less or equal to zero, then the machine precision is computed
%              by MISQP and subsequently multiplied by 1.0D+4.
%    ACCQP :   The tolerance is needed for the QP solver to perform several
%              tests, for example whether optimality conditions are satisfied
%              or whether a number is considered as zero or not. If ACCQP is
%              less or equal to zero, then the machine precision is computed
%              by MISQP and subsequently multiplied by 1.0D+4.
%    MAXIT :   Maximum number of iterations, where one iteration corresponds to
%              one evaluation of a set of gradients (e.g. 100).
%    MAXPEN :  Maximum number of successive increments of the penalty parameter 
%              without success (e.g. 50).
%    MAXUND :  Maximum number of successive iterations without improvements
%              of the iterate X (e.g. 10).
%    RESOPT :  If set to TRUE and if integer variables exist, an additional 
%              restart will be performed to check whether an improvement 
%              is possible (recommended for functions with curved narrow valleys).
%    NONMON :  Maximum number of successive iterations, which are to be 
%              considered for the non-monotone trust region algorithm.
%    IPRINT :  Specification of the desired output level.
%           0 :  No output of the program.
%           1 :  Only a final convergence analysis is given.
%           2 :  One line of intermediate results is printed in each iteration.
%           3 :  More detailed information is printed for each iteration.
%           4 :  In addition, some messages of the QP solver are displayed.
%    MODE :   The parameter allows to change some default tolerances.
%           0 :   Use default values.
%           1 :   User-provided penalty parameter SIGMA, scaling constant DELTA, 
%                 and initial trust region radii ITRC and ITRI for continuous
%                 and integer variables in RW(1), RW(2), RW(3), and RW(4).
%                 Default values are
%                   RW(1) = 10.0
%                   RW(2) = 0.01 
%                   RW(3) = 10.0
%                   RW(4) = 10.0
%    IOUT :    Integer indicating the desired output unit number, i.e., all
%              write-statements start with 'WRITE(IOUT,... '.
%    IFAIL :   The parameter shows the reason for terminating a solution
%              process. Initially IFAIL must be set to zero. On return IFAIL
%              could contain the following values:
%          -2 :   Compute new gradient values in DF and DG, see above.
%          -1 :   Compute new function values in F and G, see above.
%           0 :   Optimality conditions satisfied.
%           1 :   Termination after MAXIT iterations.
%           2 :   More than MAXUND iterations for fixing trust region.
%           3 :   More than MAXPEN updates of penalty parameter.
%           4 :   Termination at infeasible iterate.
%           5 :   Termination with zero trust region for integer variables.
%           6 :   Length of a working array is too short.
%           7 :   False dimensions, e.g., M$>$MMAX or N$>$=NMAX.
%         >90 :   QP solver terminated with an error message IFQL, 
%                 IFAIL = IFQL + 100.
%    RW(LRW),LRW :  Real working array of length LRW, and LRW must be
%                 at least 21*N + 10*M + 5*ME + 2*NONMON + 73.
%    IW(LIW),LIW :  Integer working array of length LIW, where LIW must be 
%                 at least (NINT+NBIN) + 14.
%    LW(LLW),LLW :  Logical working array of length LLW, where LLW must be
%                 at least 15.
% 
% 
%    FUNCTION AND GRADIENT EVALUATION:
% 
%    The user has to provide functions and gradients in the same program, which
%    executes also misqp, according to the following rules:
% 
%    1) Choose starting values for the variables to be optimized, and store
%       them in X.
% 
%    2) Compute objective and all constraint function values values at X and
%       store them in F and G, respectively. 
% 
%    3) Compute gradients of objective function and all constraints, and
%       store them in DF and DG, respectively. The j-th row of DG contains
%       the gradient of the j-th constraint, j=1,...,m. 
% 
%    4) Set IFAIL=0 and execute MISQP.
% 
%    5) If MISQP terminates with IFAIL=0, the internal stopping criteria are 
%       satisfied. 
% 
%    6) In case of IFAIL>0, an error occurred.
% 
%    7) If MISQP returns with IFAIL=-1, compute objective function values and
%       constraint values for all variables found in X, store them in F and G,
%       and call MISQP again. 
% 
%    8) If MISQP terminates with IFAIL=-2, compute gradient values subject to
%       variables stored in X, and store them in DF and DG. Only partial 
%       derivatives subject to the continuous variables need to be provided. 
%       Then call MISQP again.
% 
% 
% 
% 
%    AUTHOR:     O. Exler
%    -------     Department of Computer Science
%                University of Bayreuth
%                95440 Bayreuth
%                Germany
% 
% 
% 
%    VERSION:    2.2.1 (06/2008) 
%    --------
% 
% 
% ***********************************************************************
%
%m                                     integer
%me                                    integer
%n                                     integer
%nint                                  integer
%nbin                                  integer
%x (n+1)                               real
%f                                     real
%g (m+me+1)                            real
%df (n+1)                              real
%dg (m+me+1,n+1)                       real
%u (m+n+n)                             real
%xl (n+1)                              real
%xu (n+1)                              real
%b (n+1,n+1)                           real
%acc                                   real
%accqp                                 real
%maxit                                 integer
%maxpen                                integer
%maxund                                integer
%resopt                                logical
%nonmon                                integer
%iprint                                integer
%mode                                  integer
%ifail                                 integer
%rw (lrw)                              real
%lrw                                   integer
%iw (liw)                              integer
%liw                                   integer
%lw (llw)                              logical
%llw                                   integer
%
%[x,f,g,df,dg,u,xl,xu,b,ifail,rw,iw,lw] = misqp(m,me,n,nint,nbin,x,f,g,df,dg,...
%u,xl,xu,b,acc,accqp,maxit,maxpen,maxund,resopt,nonmon,iprint,mode,ifail,rw,...
%lrw,iw,liw,lw,llw)
%
%

if ifail==0
% WOKSPACE allocation
%  RW(LRW),LRW :  Real working array of length LRW, and LRW must be
%                 at least 22*N + 11*M + 
%                 + 6*ME + 2*NONMON + 85.
%  IW(LIW),LIW :  Integer working array of length LIW, where LIW must be 
%                 at least (NINT+NBIN) + 16.
%  LW(LLW),LLW :  Logical working array of length LLW, where LLW must be
%                 at least 15.

    lrw=22*n+11*m+6*me+2*nonmon+85;
    rw=zeros(1,lrw);
    liw=nint+nbin+16;
    iw=zeros(1,liw);
    llw=15;
    lw=zeros(1,llw);
end
%
%Call the MEX function
%
 [x,f,g,df,dg,u,xl,xu,b,ifail,rw,iw,lw] = misqpg(m,me,n,...
 nint,nbin,x,f,g,df,dg,u,xl,xu,b,acc,accqp,maxit,maxpen,...
 maxund,resopt,nonmon,iprint,mode,ifail,rw,lrw,iw,liw,lw,...
 llw);

if (size(x,2)==1)
    x = x';
end
if (size(g,2)==1)
    g = g';
end
if (size(df,2)==1)
    df = df';
end
% if (size(dg,2)==m+me+1)
%     dg = dg';
% end
if (size(u,2)==1)
    u = u';
end
if (size(xl,2)==1)
    xl = xl';
end
if (size(xu,2)==1)
    xu = xu';
end
if (size(rw,2)==1)
    rw = rw';
end
if (size(iw,2)==1)
    iw = iw';
end
if (size(lw,2)==1)
    lw = lw';
end


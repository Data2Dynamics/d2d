% ***********************************************************************
% 
% 
%         AN IMPLEMENTATION OF A TRUST REGION METHOD FOR SOLVING
%             MIXED-INTEGER NONLINEAR OPTIMIZATION PROBLEMS 
% 
% 
%    MISQP SOLVES THE MIXED-INTEGER NONLINEAR PROGRAM (MINLP)
% 
%              MINIMIZE        F(X,Y)
%              SUBJECT TO      G(J)(X,Y)   =  0  , J=1,...,ME
%                              G(J)(X,Y)  >=  0  , J=ME+1,...,M
%                              XL  <=  X  <=  XU
%                              YL  <=  Y  <=  YU
% 
%    WHERE X IS A REAL AND Y AN INTEGER VARIABLE VECTOR.
% 
%    THE FORTRAN SUBROUTINE IS AN IMPLEMENTATION OF A MODIFIED SEQUENTIAL 
%    QUADRATIC PROGRAMMING (SQP) METHOD. UNDER THE ASSUMPTION THAT INTEGER 
%    VARIABLES HAVE A 'SMOOTH' INFLUENCE ON THE MODEL FUNCTIONS, I.E., THAT 
%    FUNCTION VALUES DO NOT CHANGE TOO DRASTICALLY WHEN IN- OR DECREMENTING 
%    AN INTEGER VALUE, SUCCESSIVE QUADRATIC APPROXIMATIONS ARE APPLIED.
%    THE ALGORITHM IS STABILIZED BY A TRUST REGION METHOD WITH YUAN'S SECOND 
%    ORDER CORRECTIONS. 
% 
%    IT IS NOT ASSUMED THAT THE MIXED-INTEGER PROGRAM IS RELAXABLE. IN OTHER 
%    WORDS, FUNCTION VALUES ARE REQUIRED ONLY AT INTEGER POINTS. THE HESSIAN 
%    OF THE LAGRANGIAN FUNCTION IS APPROXIMATED BY BFGS UPDATES SUBJECT TO 
%    THE CONTINUOUS VARIABLES. DIAGONAL SECOND ORDER INFORMATION SUBJECT TO 
%    THE INTEGER VARIABLES IS OBTAINED BY A DIFFERENCE FORMULA.
% 
% 
%    USAGE:
% 
%       [x,f,g,df,dg,u,xl,xu,b,ifail,rw,iw,lw] = mi_sqp(m,me,mme,n,nint,nbin,x,f,g,...
% 	df,dg,u,xl,xu,b,acc,accqp,maxit,maxpen,maxund,iprint,mode,ifail,rw,lrw,iw,...
% 	liw,lw,llw)
% 
% 
%    ARGUMENTS: (NMAX=N+1, MMAX=M+ME+1)
% 
%    M :       TOTAL NUMBER OF CONSTRAINTS.
%    ME :      NUMBER OF EQUALITY CONSTRAINTS.
%    N :       NUMBER OF OPTIMIZATION VARIABLES, CONTINUOUS AND INTEGER ONES.
%    NINT :    NUMBER OF INTEGER VARIABLES, MUST BE LESS THAN OR EQUAL TO N.
%    NBIN :    NUMBER OF BINARY VARIABLES, MUST BE LESS THAN OR EQUAL TO N.
%    MME :     MUST BE EQUAL TO M+ME
%    X(N+1) :  INITIALLY, X HAS TO CONTAIN STARTING VALUES. ON RETURN, X IS
%              REPLACED BY THE CURRENT ITERATE. IN THE DRIVING PROGRAM
%              THE ROW DIMENSION OF X HAS TO BE EQUAL TO NMAX.
%    F :       F CONTAINS THE ACTUAL OBJECTIVE FUNCTION VALUE EVALUATED AT X.
%    G(MME+1) : G CONTAINS THE ACTUAL CONSTRAINT FUNCTION VALUES AT X.
%    DF(N+1) : DF CONTAINS THE GRADIENT OF THE OBJECTIVE FUNCTION.
%    DG(MME+1,N+1) : DG CONTAINS GRADIENTS OF CONSTRAINTS IN THE FIRST M ROWS. 
%              IN THE DRIVING PROGRAM THE ROW DIMENSION OF DG MUST BE EQUAL TO 
%              MMAX.
%    U(M+N+N) :  ON RETURN, U CONTAINS MULTIPLIERS. THE FIRST M LOCATIONS CONTAIN
%              THE MULTIPLIERS OF THE M NONLINEAR CONSTRAINTS, THE SUBSEQUENT
%              N LOCATIONS THE MULTIPLIERS SUBJECT TO THE LOWER BOUNDS, AND THE
%              FINAL N LOCATIONS THE MULTIPLIERS SUBJECT TO THE UPPER BOUNDS.
%              AT AN OPTIMAL SOLUTION, ALL MULTIPLIERS WITH RESPECT TO
%              INEQUALITY CONSTRAINTS SHOULD BE NONNEGATIVE.
%    XL(N+1),XU(N+1) : ON INPUT, THE ONE-DIMENSIONAL ARRAYS XL AND XU MUST
%              CONTAIN THE UPPER AND LOWER BOUNDS OF THE VARIABLES.
%    B(N+1,N+1) : ON RETURN, B CONTAINS THE LAST COMPUTED APPROXIMATION
%              OF THE HESSIAN MATRIX OF THE LAGRANGIAN FUNCTION. 
%              IN THE DRIVING PROGRAM, THE ROW DIMENSION OF C HAS TO BE EQUAL
%              TO NMAX.
%    ACC :     THE USER HAS TO SPECIFY THE DESIRED FINAL ACCURACY
%              (E.G. 1.0D-7). THE TERMINATION ACCURACY SHOULD NOT BE SMALLER
%              THAN THE ACCURACY BY WHICH GRADIENTS ARE COMPUTED. IF ACC IS
%              LESS OR EQUAL TO ZERO, THEN THE MACHINE PRECISION IS COMPUTED 
%              BY MISQP AND SUBSEQUENTLY MULTIPLIED BY 1.0D+4.
%    ACCQP :   THE TOLERANCE IS NEEDED FOR THE QP SOLVER TO PERFORM SEVERAL 
%              TESTS, FOR EXAMPLE WHETHER OPTIMALITY CONDITIONS ARE SATISFIED
%              OR WHETHER A NUMBER IS CONSIDERED AS ZERO OR NOT. IF ACCQP IS
%              LESS OR EQUAL TO ZERO, THEN THE MACHINE PRECISION IS COMPUTED 
%              BY MISQP AND SUBSEQUENTLY MULTIPLIED BY 1.0D+4.
%    MAXIT :   MAXIMUM NUMBER OF ITERATIONS, WHERE ONE ITERATION CORRESPONDS TO 
%              ONE EVALUATION OF A SET OF GRADIENTS (E.G. 100).
%    MAXPEN :  MAXIMUM NUMBER OF SUCCESSIVE INCREMENTS OF THE PENALTY PARAMETER
%              WITHOUT SUCCESS (E.G. 50).
%    MAXUND :  MAXIMUM NUMBER OF SUCCESSIVE ITERATIONS WITHOUT IMPROVEMENTS OF
%              THE ITERATE X (E.G. 10).
%    IPRINT :  SPECIFICATION OF THE DESIRED OUTPUT LEVEL.
%       IPRINT = 0 :  NO OUTPUT OF THE PROGRAM.
%       IPRINT = 1 :  ONLY A FINAL CONVERGENCE ANALYSIS IS GIVEN.
%    MODE :    THE PARAMETER SPECIFIES THE DESIRED VERSION OF MISQP.
%              FOR MODE<4, DERIVATIVE APPROXIMATIONS SUBJECT TO THE INTEGER
%              VARIABLES ARE CALCULATED BY MISQP, AND ARE OTHERWISE TO BE
%              PROVIDED BY THE USER. 
%       MODE = 0 :    NORMAL EXECUTION (REVERSE COMMUNICATION!).
%       MODE = 1 :    THE USER PROVIDES AN INITIAL GUESS FOR THE PENALTY
%                     PARAMETER SIGMA, THE SCALING CONSTANT DELTA, AND THE 
%                     INITIAL TRUST REGION RADII ITRC AND ITRI FOR CONTINUOUS 
%                     AND INTEGER VARIABLES IN RW(1), RW(2), RW(3), AND RW(4). 
%                     DEFAULT VALUES ARE     
%                         RW(1) = 1.0D+1
%                         RW(2) = 1.0D-2
%                         RW(3) = 1.0D+1      
%                         RW(4) = 1.0D+1     
%    IFAIL :   THE PARAMETER SHOWS THE REASON FOR TERMINATING A SOLUTION
%              PROCESS. INITIALLY IFAIL MUST BE SET TO ZERO. ON RETURN IFAIL
%              COULD CONTAIN THE FOLLOWING VALUES:
%       IFAIL =-2 :   COMPUTE GRADIENT VALUES SUBJECT TO THE VARIABLES STORED IN
%                     X, AND STORE THEM IN DF AND DG. THEN CALL MISQP AGAIN, 
%                     SEE BELOW.
%       IFAIL =-1 :   COMPUTE OBJECTIVE FUNCTION AND ALL CONSTRAINT VALUES SUBJECT
%                     THE VARIABLES FOUND IN X, AND STORE THEM IN F AND G. 
%                     THEN CALL MISQP AGAIN, SEE BELOW.
%       IFAIL = 0 :   THE OPTIMALITY CONDITIONS ARE SATISFIED.
%       IFAIL = 1 :   THE ALGORITHM HAS BEEN STOPPED AFTER MAXIT ITERATIONS.
%       IFAIL = 2 :   THERE ARE TOO MANY SUBITERATIONS TO FIX THE TRUST REGION,
%                     MORE THAN MAXUND. 
%       IFAIL = 3 :   THERE ARE MORE THAN MAXPEN CORRECTIONS OF THE PENALTY 
%                     PARAMETER WITHOUT IMPROVEMENTS.
%       IFAIL = 4 :   SEARCH ALGORITHM STOPPED, BUT THE FINAL ITERATE IS STILL 
%                     INFEASIBLE.
%       IFAIL = 5 :   SEARCH ALGORITHM STOPPED AND THE TRUST REGION RADIUS OF THE 
%                     INTEGER VARIABLES IS 0. 
%       IFAIL = 6 :   LENGTH OF A WORKING ARRAY IS TOO SHORT.
%       IFAIL = 7 :   THERE ARE FALSE DIMENSIONS, FOR EXAMPLE M>MMAX OR N>=NMAX.
%       IFAIL > 90 :  THE SOLUTION OF THE QUADRATIC PROGRAMMING SUBPROBLEM HAS BEEN 
%                     TERMINATED WITH AN ERROR MESSAGE AND IFAIL IS SET TO IFQL+100.
%       RW(LRW),LRW : RW IS A REAL WORKING ARRAY OF LENGTH LRW. LRW MUST BE 
%                     AT LEAST lrw=100+19*max(1,m)+23*max(1,n)+max(1,m)*max(1,n).
%       IW(LIW),LIW : THE USER HAS TO PROVIDE WORKING SPACE FOR AN INTEGER ARRAY
%                     OF LENGTH LIW. LIW MUST BE AT LEAST max(nint+nbin,1)+15.
%       LW(LLW),LLW : THE USER HAS TO PROVIDE WORKING SPACE FOR A LOGICAL ARRAY
%                     OF LENGTH LLW. LLW MUST BE AT LEAST 15.
% 
% working arrays
% 
% lrw=100+19*max(1,m)+23*max(1,n)+max(1,m)*max(1,n);
% rw=zeros(1,lrw);
% liw=max(nint+nbin,1)+15;
% iw=zeros(1,liw);
% llw=15;
% lw=zeros(1,llw);
% 
% 
%    FUNCTION AND GRADIENT EVALUATION:
% 
%    THE USER HAS TO PROVIDE FUNCTIONS AND GRADIENTS IN THE SAME PROGRAM, WHICH
%    EXECUTES ALSO MISQP, ACCORDING TO THE FOLLOWING RULES:
% 
%    1) CHOOSE STARTING VALUES FOR THE VARIABLES TO BE OPTIMIZED, AND STORE
%       THEM IN X.
% 
%    2) COMPUTE OBJECTIVE AND ALL CONSTRAINT FUNCTION VALUES VALUES AT X AND
%       STORE THEM IN F AND G, RESPECTIVELY. 
% 
%    3) COMPUTE GRADIENTS OF OBJECTIVE FUNCTION AND ALL CONSTRAINTS, AND
%       STORE THEM IN DF AND DG, RESPECTIVELY. THE J-TH ROW OF DG CONTAINS
%       THE GRADIENT OF THE J-TH CONSTRAINT, J=1,...,M. IN CASE OF MODE<4, 
%       ONLY PARTIAL DERIVATIVES SUBJECT TO THE CONTINUOUS VARIABLES NEED TO 
%       BE PROVIDED.
% 
%    4) SET IFAIL=0 AND EXECUTE MISQP.
% 
%    5) IF MISQP TERMINATES WITH IFAIL=0, THE INTERNAL STOPPING CRITERIA ARE 
%       SATISFIED. 
% 
%    6) IN CASE OF IFAIL>0, AN ERROR OCCURRED.
% 
%    7) IF MISQP RETURNS WITH IFAIL=-1, COMPUTE OBJECTIVE FUNCTION VALUES AND
%       CONSTRAINT VALUES FOR ALL VARIABLES FOUND IN X, STORE THEM IN F AND G,
%       AND CALL MISQP AGAIN. 
% 
%    8) IF MISQP TERMINATES WITH IFAIL=-2, COMPUTE GRADIENT VALUES SUBJECT TO
%       VARIABLES STORED IN X, AND STORE THEM IN DF AND DG. IN CASE OF MODE<4, 
%       ONLY PARTIAL DERIVATIVES SUBJECT TO THE CONTINUOUS VARIABLES NEED TO 
%       BE PROVIDED. THEN CALL MISQP AGAIN.
% 
% 
% ***********************************************************************
%             [NAGWare Gateway Generator]
%
%Copyright (c) 1993-97 by the Numerical Algorithms Group Ltd 2.0a
%
%
% maximal n=100 and m=100 !!!!!!!!!!!
%
%m                                     integer
%me                                    integer
%mme                                   integer
%n                                     integer
%nint                                  integer
%nbin                                  integer
%x (n+1)                               real
%f                                     real
%g (mme+1)                             real
%df (n+1)                              real
%dg (mme+1,n+1)                        real
%u (m+n+n)                             real
%xl (n+1)                              real
%xu (n+1)                              real
%b (n+1,n+1)                           real
%acc                                   real
%accqp                                 real
%maxit                                 integer
%maxpen                                integer
%maxund                                integer
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
%[x,f,g,df,dg,u,xl,xu,b,ifail,rw,iw,lw] = mi_sqp(m,me,mme,n,nint,nbin,x,f,g,...
%df,dg,u,xl,xu,b,acc,accqp,maxit,maxpen,maxund,iprint,mode,ifail,rw,lrw,iw,...
%liw,lw,llw)
%
%
 function [x,f,g,df,dg,u,xl,xu,b,ifail,rw,iw,lw] = mi_sqp(m,me,mme,n,nint,...
 nbin,x,f,g,df,dg,u,xl,xu,b,acc,accqp,maxit,maxpen,maxund,iprint,mode,ifail,...
 rw,lrw,iw,liw,lw,llw)
%
%
%
%Call the MEX function
%
 [x,f,g,df,dg,u,xl,xu,b,ifail,rw,iw,lw] = mi_sqpg(m,me,...
 mme,n,nint,nbin,x,f,g,df,dg,u,xl,xu,b,acc,accqp,maxit,...
 maxpen,maxund,iprint,mode,ifail,rw,lrw,iw,liw,lw,llw);
%

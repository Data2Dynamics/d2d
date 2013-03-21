function[sol,ierr,output,history,grad,diagnostic]= ...
                                    STRSCNE(x,F,tol,l,u,parms,Fjac)
%  
%   Globally convergent  solver for constrained nonlinear
%   systems of equations
%
%           F(x)=0         l<=x<=u
%
%   where  F: R^n --> R^n. 
%   The algorithm combines Newton method and an elliptical trust-region 
%   procedure. The  merit function  used is f(x)=0.5*norm(F(x))^2.
%   The elliptical trust-region is defined employing a scaling 
%   diagonal matrix D and the trust-region subproblem is approximately 
%   solved by the dogleg method.  Strictly feasible iterates are generated.
%   For more details, see
%
%   S.Bellavia, M.Macconi, B.Morini,
%       `` An affine scaling trust region method approach to  
%       bound-constrained  nonlinear systems '', 
%       Applied Numerical Mathematics, to appear. 
% 
%   S.Bellavia, M.Macconi, B.Morini, 
%      STRSCNE: A Scaled Trust Region Solver for
%      Constrained Nonlinear Equations
%      Technical Report, 2002
%
%   email: bellavia@ciro.de.unifi.it
%          morini@ciro.de.unifi.it
%
%      Requires  dmatrici.m, diffjac.m, dogleg.m
%
%
%
%  function [sol,ierr,output]=STRSCNE(x,F,tol,l,u,parms)
%  -----------------------------------------------------
%  
%  The Jacobian matrix of F is approximated using 
%  finite-differences.
%
%  inputs:
%  
%     x     = initial iterate x_0.
%     F     = function that accepts input vector X and returns the 
%             vector F(X).
%     tol   = [atol,rtol] absolute/relative error tolerances used 
%             in the stopping criterion:
%                    norm(F(x))<=atol+rtol*norm(F(x_0)).
%     l,u   = vectors containing the lower and upper bounds on the
%             variables.
%             set l(i) = -Inf if x(i) is unbounded below; 
%             set u(i) = Inf if x(i) is unbounded above.
%     parms = [maxit,maxnf,delta,outflag]
%             maxit=  maximum number of nonlinear iterations;
%             maxnf=  maximum number of F-evaluations (F-evaluations 
%                     due to the Jacobian approximation are not included);
%             delta=  choice of the initial trust region radius Delta. 
%                     -1  then Delta=1;
%                     -2  then Delta=norm( inv(D(x_0))grad(f(x_0)) );
%                     >0  then Delta=delta.
%             outflag= printlevel 
%                     0   The final summary output is printed:
%                         number of performed iterations,  
%                         norm(F(sol))
%                         number of performed F-evaluations (F-evaluations 
%                         due to the Jacobian approximation are not 
%                         included);
%                     >0  one line of summary output for each iteration 
%                         is printed:
%                         iteration number k;
%                         norm(F(x_k));
%                         number of  trust region radius reductions;  
%                         step back taken to stay feasible;
%                         ratio of two successive nonlinear residuals; 
%                         direction used: 
%                                       nw: Newton direction
%                                       c:  scaled gradient direction 
%                                       d:  dogleg direction (when 
%                                           different from nw) 
%  output:
%     sol    =   final estimate of the solution.
%     ierr  
%         =  0   upon successful termination;
%         =  1   the limiting number of iterations has been reached; 
%         =  2   the limiting number of F-evaluations has been reached;
%         =  3   the trust region radius Delta has become too small 
%                (Delta<sqrt(eps));
%         =  4   no improvement for the nonlinear reasidual could be obtained:
%                abs(norm(F(x_k))-norm(F(x_{k-1})))<=100*eps*norm(F(x_k));
%         =  5   the sequence has approached a minumum of f the box:
%                norm(inv(D(x_k)*grad(f(x_k)))<100*eps. 
%         =  6   an overflow would be generated when computing the scaling 
%                matrix D since the sequence is approaching a bound.       
%     output =   vector containing:
%                number of performed iterations;
%                number of performed F-evaluations (F-evaluations due to 
%                the Jacobian approximation are not included); 
%                norm(F(sol));
%                norm( inv(D(sol))grad(f(sol)) );  
%                total number of reductions of the trust region radius. 
%
%
%
%
%
%  function [sol,ierr,output]=STRSCNE(x,F,tol,l,u,parms,Fjac)
%  ---------------------------------------------------------- 
%
%  Solves as above with the Jacobian matrix of F evaluated analitically 
%  by the user-supplied function Fjac. Function Fjac(X) must  
%  returns the Jacobian matrix of the function F evaluated at X. 
%
%  function [sol,ierr,output,history]=STRSCNE(x,F,tol,l,u,parms,..) 
%  
%  Returns the  matrix history that describes the convergence 
%  history of STRSCNE:
%  Each row of history contains:
%                        iteration number k; 
%                        norm(F(x_k));
%                        number of  trust region radius reductions;  
%                        step back taken to stay feasible;
%                        ratio of two successive nonlinear residuals. 
% 
%
%
%
%  function [sol,ierr,output,history,grad]=STRSCNE(x,F,tol,l,u,parms,..)
%  ---------------------------------------------------------------------
%
%  Returns the gradient grad of the merit function f at sol.
%
%
%
%
%  function [sol,ierr,output,history, grad,diagnostic]
%                                     =STRSCNE(x,F,tol,l,u,parms,..)
%  ----------------------------------------------------------------- 
%
%  Returns some diagnostic information: 
%                 diagnostic(1)=rank of the Jacobian matrix at sol, 
%                               computed by the Matlab function rank;
%                 diagnostic(2:n+1) singular values of the Jacobian 
%                               matrix at sol, computed by the Matlab 
%                               function svd.
%
%
% Internal parameters:
%   r=0.1, t=0.25        : used for accuracy requirements;
%   thetal=0.99995       : used to ensure strictly feasible  iterates; 
%   delta1=0.25,delta2=2 : used to update the trust-region size; 
%   w=0.75               : parameter that governs the increase of the 
%                          trust-region size.
%
%   Initialization
n=length(x);  
ierr=0; nridut=0;
itc=0;  %iterations counter
nvf=0;  %counter of the F-evaluations 
%
for i=1:n
    if (x(i)<=l(i) | x(i)>=u(i))
%
%   The given initial guess is not feasible,
%
        disp('ATTENTION: the initial guess is not feasible')
          return
    end
    
end
fx=feval(F,x); nvf=nvf+1; fnrm=norm(fx); fnrm2=fnrm^2;
%
%   Parameter settings
%
epsilon=100*eps; 
Deltamin=sqrt(eps);       
atol=tol(1); rtol=tol(2); stoptol=atol+rtol*fnrm;
maxit=parms(1); maxnf=parms(2); Delta=parms(3); outflag=parms(4);
%   
%   Internal parameters
%
r=0.1; t=0.25; w=0.75; delta1=0.25; delta2=2; thetal=0.99995;  
%
%   Initial output
if nargout>=4
     history(1,:)=[itc,fnrm, 0, 0 ,0 ];
end
if outflag>0
    disp(sprintf('\n'))
    disp(sprintf('%s %s %s %s %s %s %s %s %s %s %s %s', ...
    blanks(6),'it',blanks(4),'||F||_2',blanks(2),'rid_step', blanks(2),'alpha', ...
    blanks(7),'ratio', blanks(4),'direc'))
end 
%
%   Iteration
%
while (fnrm>stoptol & itc< maxit & nvf<maxnf)
  itc=itc+1; fnrm0=fnrm;
%
%   Jacobian evaluation
%
     if nargin==6
       jac= diffjac(x,F,fx,l,u);
     else
       jac=feval(Fjac,x);
     end    
     grad= jac'*fx;
%
%   Calculation of the scaling matrices D (d), inv(D) (dm1),
%   and inv(D)^2 (dm2) and of the matrix jjac=jac'*jac
 
    [d,dm1,dm2,ierr]=dmatrici(x,grad,l,u); 
    if ierr==6
       break;
    end
    jjac=jac'*jac;
%  
  
    dm1grad=dm1.*grad; ndm1grad=norm(dm1grad);
    if ndm1grad<epsilon
       ierr=5;
       break;
    end     
%
%   Initial trust region radius
%
    if itc==1
       if Delta==-1
          Delta=1; 
       else  if Delta==-2
                Delta=norm(dm1.*grad);
             end
       end
    end      
%
%   Newton step
%
    sn=jac\(-fx);
%
%   Computation of the minimizer (pc) of the quadratic model 
%   along dm2*grad 
%
    dm2grad=dm2.*grad; 
    vert=(ndm1grad/norm(jac*dm2grad))^2;
%
%   
  pc=-vert*dm2grad;
  pcv=pc;  npc=norm(pc);
%  
%  computation of D*sn (snc) and D*pc (pcc)  
%
  snc=d.*sn; pcc=d.*pc;
%
%  trust-region strategy
%
  rhof=0; nridu=-1;
  while ( rhof<t & Delta>Deltamin)    
         nridu=nridu+1;           
%
%       Computation of the Cauchy point
%
         if(norm(pcc)>Delta)
            pcv=-Delta*dm2grad/ndm1grad;
         end
%
%        Computation of the truncated Cauchy point
%
         pciv=pcv;
         npcv=norm(pcv);
         for i=1:n
            if pcv(i)~=0 
               alp(i)=max((l(i)-x(i))/pcv(i),(u(i)-x(i))/pcv(i));
            else
               alp(i)=Inf;
            end
         end        
         alpha=min(alp);
         if (alpha<=1)
             pciv=max(thetal,1-npcv)*alpha*pcv;
         end
         if(norm(snc)<=Delta)
%
%         The Newton step is the solution of the trust-region subproblem.
%         Computation of the truncated Newton step. 
%
            p=sn;
            np=norm(p);
            for i=1:n
              if p(i)~=0 
                alp(i)=max((l(i)-x(i))/p(i),(u(i)-x(i))/p(i));
              else
                alp(i)=Inf;
              end
            end        
            alpha=min(alp);
            if (alpha<=1)
                p=max(thetal,1-np)*alpha*p;
            end
            step='nw'; alpha=norm(p)/np;
         else
            if (norm(pcc)>=Delta)
%
%         The Cauchy step is the solution of the trust-region subproblem
% 
               p=pciv;
               step='c'; alpha=norm(pciv)/npcv;
            else
%
%          Dogleg method
%                           
               dc=dogleg(snc,pcc,Delta);
               pd=dm1.*dc;
               npd=norm(pd); p=pd;
               for i=1:n
                 if pd(i)~=0 
                    alp(i)=max((l(i)-x(i))/pd(i),(u(i)-x(i))/pd(i));
                 else
                    alp(i)=Inf;
                 end
               end        
               alpha=min(alp);
               if (alpha<=1)
                   p=max(thetal,1-npd)*alpha*pd;
               end
               step='d'; alpha=norm(p)/npd;
            end
         end
%
%         Accuracy requirements
%
         rhoc=(grad'*p+0.5*p'*jjac*p)/ ...
              (grad'*pciv+0.5*pciv'*jjac*pciv);
         if (rhoc<r)    
%
%          switch to the truncated Cauchy step
%         
                p=pciv;
                step='c'; alpha=norm(pciv)/npcv;
         end
         xpp=x+p';
         fxpp=feval(F,xpp);
         nvf=nvf+1;
         fnrmxpp=norm(fxpp);
         rhof=(fnrmxpp^2-fnrm2)*0.5/ ...
              (grad'*p+0.5*p'*jjac*p);
         Deltas=Delta;
         Delta=min(delta1*Delta, 0.5*norm(d.*p));  
    end
%
if (Delta <= Deltamin & rhof<t)
      nridu=nridu+1;
      ierr=3;
      break
  end 
  Delta=Deltas;
%  
%     Updating of the iterate.
%
  x=xpp;
  fx=fxpp; 
  fnrm=fnrmxpp; fnrm2=fnrm^2; rat=fnrm/fnrm0; nridut=nridut+nridu;
%
%  Storing and printing the iteration's summary.
% 
  if nargout>-4
     history(itc+1,:)=[itc, fnrm, nridu, alpha, rat];
  end
  if outflag>0  
	disp(sprintf('%s %6.0f %s %10.5e %s %3.0f %s %10.5e %s %10.5e  %s %s', ...
	blanks(2),  itc, blanks(2),fnrm, blanks(2),nridu, blanks(2), ...
	alpha, blanks(2),rat,blanks(2), step))
  end
  if (abs(fnrm-fnrm0)<=epsilon*fnrm & fnrm>stoptol)  
      ierr=4;
      break
      return
  end
%
%  Updating of the trust-region size
%  
  if (rhof>w & rhoc>r)   
      Delta=max(Delta, delta2*norm(d.*p));
  end
end
%
%  Final output
% 
sol=x;
if nargin==6 
  jac= diffjac(x,F,fx,l,u);
else
  jac=feval(Fjac,x);
end    
grad= jac'*fx;
dkm1=ones(n,1);
for i=1:n
    if (grad(i)<0) 
        if (u(i)~=Inf)
             diff=u(i)-x(i);
             dkm1(i)=sqrt(diff);
        end  
    else
        if (l(i)~=-Inf)
             diff=x(i)-l(i);
             dkm1(i)=sqrt(diff);
        end 
    end
end
ndm1grad=norm(dkm1.*grad);    
output=[itc,nvf,fnrm,ndm1grad,nridut];
if nargout==6
   diagnostic(1)=rank(jac);
   diagnostic(2:n+1)=svd(jac);
end   
    

if (ierr==0 & fnrm>stoptol)
    if itc==maxit
            ierr=1;
     else       
           ierr=2;
     end
end  
if (ierr>=1)
     disp(sprintf('%s %1.0f','FAILURE, ierr= ',ierr))
     return     
end 


disp(sprintf('\n'))
disp(sprintf('%s %6.0f','iterations = ', itc)) 
disp(sprintf('%s %10.5e','||F||_2= ',fnrm))
disp(sprintf('%s %6.0f','F-evaluations (no Jacobian) = ', nvf)) 
   


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function[jac]=diffjac(x,f,f0,l,u);
%
%   compute a the Jacobian matrix of f at x by forward finite difference.s
%
%    function[jac]=diffjac(x,f,f0,l,u);
% 
%  inputs:
%  x, f   = point and function.
%  f0      = f(x), preevaluated.
%  l,u     = constraints.
%
%  output 
%     jac  = approximated Jacobian matrix. 
%
n=length(x);
epsnew=sqrt(eps);
for j=1:n
%    
%      choice of the steplenght for the forward differences
% 
   if  x(j)==0
       h=epsnew;
    else
       h=epsnew*sign(x(j))*max(abs(x(j)),norm(x,1)/n);
   end
   xhj=x(j)+h; 
%
%
   if (xhj<l(j) |xhj >u(j))
%
%       the new point xhj is not feasible. In this case the backward 
%       difference is used
%
       h=-h;
       xhj=x(j)+h;
       if (xhj<l(j) |xhj >u(j))
           disp('Function diffjac: loosing of feasibility using both backward and forward differences')
           stop
       end
   end
   xh=x;  xh(j)=xhj;
   f1=feval(f,xh);
       jac(:,j)=(f1-f0)/h;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[dk,dkm1,dkm2,ind]=dmatrici(x,grad,l,u);

%
%  
%   Compute the scaling matrix D and the metrices inv(D), (D)^(-2)
%
%      function[dk,dkm1,dkm2,ind]=dmatrici(x,grad,l,u);
%
% input 
% x  = point.
% grad = gradient of f, preevaluated.
% l,u     = constraints.
% 

% output 
% dk   =  array such that dk=diag(D(x)).
% dkm1 =  array such that dkm1=diag(inv(D(x))).
% dkm2 =  array such that dkm1=diag((D(x))^-2).
% ind  =  0 upon successful termination.
%      =  6 an overflow would be generated when computing 
%         the scaling matrix D.
%
ind=0;
n=length(x);
dk=ones(n,1); dkm1=ones(n,1);  dkm2=ones(n,1); 
for i=1:n
    if (grad(i)<0) 
        if (u(i)~=Inf)
             diff=u(i)-x(i);
             sqdiff=sqrt(diff);
             if diff>=(1/realmax)
                  dk(i)=1/sqdiff;
                  dkm1(i)=sqdiff;
                  dkm2(i)=diff;
             else
               ind=6;
               return
            end 
        end
    else
        if (l(i)~=-Inf)
             diff=x(i)-l(i);
             sqdiff=sqrt(diff);
          if diff>=(1/realmax)
               dk(i)=1/sqdiff;
               dkm1(i)=sqdiff;
               dkm2(i)=diff;
          else
               ind=6;
               return
          end 
        end
    end
end
return
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
function[sol]=dogleg(s,p,Delta)
%
%    Compute the dogleg step.
%
%      function[sol]=dogleg(s,p,Delta)
% input:
% s     = Newton step
% p     = Cauchy point
% Delta = trust-region size
%
% output:
% sol   = dogleg step
%
pnorm2=norm(p)^2;
a=norm(s-p)^2;
b=(p'*s-pnorm2);
c=pnorm2-Delta^2;
l1=(-b+sign(-b)*sqrt(b^2-a*c))/a;
l2=c/(l1*a);
lambda=max(l1,l2);
sol=p+lambda*(s-p);
return

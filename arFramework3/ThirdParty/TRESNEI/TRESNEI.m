function [sol,ierr,output] = TRESNEI(x,e_i,fun,l,u,options,varargin)    
%  TRESNEI solves systems of nonlinear equalities and inequalities
%
%  TRESNEI implements a trust-region Gauss-Newton method for 
%  bound-constrained least-squares problem:
%                                 
%                        min   || F(x) ||_2^2
%                      l<=x<=u             
%  
%  where F:R^n -> R^m.
%
%  Function TRESNEI solves: 
%   i)   bound-constrained square and nonsquare systems of nonlinear equations
%           C_E(x)=0,         C_E:R^n-> R^me,
%           l<=x<=u, 
%
%   ii)  nonlinear least-squares 
%            min || C_E(x) ||_2^2,   C_E:R^n -> R^me,
%          l<=x<=u   
%
%   iii) systems of nonlinear equalities and inequalities
%            C_E(x)=0,
%            C_I(x)<=0,        C_E:R^n -> R^me,   C_I:R^n -> R^mi.
%            l<=x<=u,
%
%
%  The bound-constrained least-squares problem solved internally is
%       
%                                      ||    _                    _   || 2
%                                      ||   |        C_E(x)        |  ||
%     min   || F(x) ||_2^2 =     min   ||   |                      |  ||,    
%   l<=x<=u                    l<=x<=u ||   |  0.5*max(0,C_I(x))^2 |  ||
%                                      ||    _                    _   ||2
%  where F:R^n -> R^(me+mi).
%       
%  Input: 
%   
%   x       = initial iterate such that l<=x<=u. 
%   e_i     = vector [me,mi] that holds the length me of vector C_E  and 
%             the length mi of vector C_I. 
%   fun     = function that accepts input vector x and returns a vector F 
%                 function F= fun(x)
%             The vector F is defined as: 
%             i)   the  column vector residual function C_E of length me if the 
%                  problem is a nonlinear system or a nonlinear least-squares;
%                                     _       _
%             ii)  the column vector |   C_E   |  of length me+mi in case of systems 
%                                    |_  C_I  _|
%                  of nonlinear equalities and inequalities (the equalities, if any, 
%                  must be the first components of the returned vector);
%     
%             If the Jacobian of C can be computed and the Jacobian option 
%             options.jacobian is  'on', then the function  fun 
%             must return the Jacobian matrix J at x in a second output
%             argument:
%                 function [F,J]= fun(x)
%             Checking the value of nargout, the computation of J can be avoided   
%             when the function fun is called with only one output argument.
%   l,u     = vectors containing the lower and upper bounds on the variables. 
%             set l(i) = -1.e20 if x(i) is unbounded below;  
%             set u(i) =  1.e20 if x(i) is unbounded above. 
%   options = (optional)  structure with user's parameters. The special call 
%                       options= TRESNEI()
%               returns a structure with default values.  Valid fields are:
%              .max_itns : maximum number of nonlinear iterations 
%                          (default value is 1000);   
%              .max_feval: maximum number of function F-evaluations 
%                          (F-evaluations for Jacobian finite difference 
%                           approximation are not included)
%                          (default value is 1000); 
%              .delta    : initial trust region radius (default value is 1)
%              .tol_F    : termination tolerance on norm(F, inf) 
%                          (default value is 1.e-6)
%              .tol_opt  : termination tolerance on the first order optimality 
%                          condition (default value is 1.e-6)
%              .jacobian : if 'on' the user's defined Jacobian provided by fun 
%                          is used, otherwise the jacobian is approximated by finite 
%                          differences (default value is 'on')
%              .output   : print level  (default value is 1)
%                         .output=0   a final summary output is printed: 
%                                     number of performed iterations,   
%                                     number of performed F-evaluations 
%                                       (F-evaluations for Jacobian approximation 
%                                        are not included); 
%                                     norm of the final residual ||F||_2;
%                                     first-order optimality measure.
%                         .output<0   no output is displayed.
%                         .output>0   one line of summary output for each iteration is 
%                                     printed: 
%                                     iteration number k; 
%                                     ||F(x_k)||_2;
%                                     first-order optimality measure;
%                                     trust-region solution:  N (Newton direction), 
%                                                       C (scaled gradient direction),  
%                                                       D (dogleg direction) ; 
%                                     trust-region radius; 
%                                     norm of the step; 
%                                     direction of the step: TR (trust region step), 
%                                               PTR (projected trust-region step), 
%                                               CC (convex combination of the 
%                                                   trust-region step and the Cauchy step).
%   varargin= additional parameters P1,P2,... passed to the function fun as 
%             fun(x,P1,P2...).
% 
%   Output: 
%     sol   =  final estimate of the solution. 
%     ierr  =  output flag. Possible values are: 
%               0   the nonlinear residual condition on norm(F,inf) is met;
%               1   the first-order optimality condition is met;
%               2   the limiting number of iterations has been reached;  
%               3   the limiting number of F-evaluations has been reached; 
%               4   the trust region radius Delta has become too small (Delta<eps); 
%               5   no improvement for the nonlinear reasidual could be obtained: 
%                   | ||F(x_k)||_2 - ||F(x_{k-1})||_2 |<=100*eps*||F(x_k)||_2; 
%    output =  structure with fields: 
%                .itns = number of performed iterations; 
%                .feval = number of performed F-evaluations (F-evaluations due to  
%                       the Jacobian approximation are not included);  
%                .Fnorm = ||F(sol)||_2; 
%                .optimality = first-order optimality  measure;   
%                .violations = constraint violations norm(max(C_I(sol),0),inf);
%                
%   References: 
%   
%   [1] B.Morini, M.Porcelli, "TRESNEI, a Matlab trust-region solver for systems 
%       of nonlinear equalities and inequalities", Computational
%       Optimization and Applications, 51:1, pp. 27-49, 2012. 
%   [2] M.Macconi, B.Morini, M.Porcelli, "Trust-region quadratic methods for 
%       nonlinear systems of mixed equalities and inequalities", 
%       Applied Numerical Mathematics, 59, pp. 859-876, 2009.
%   [3] M.Macconi, B.Morini, M.Porcelli, "A Gauss-Newton method for solving 
%       bound-constrained underdetermined nonlinear systems", 
%       Optimization Methods and Software, 24, pp. 219-235, 2009.
% 
%   Last Update: March 23, 2012
% 
defaultOpts.max_itns  = 1000;   
defaultOpts.max_feval = 1000;
defaultOpts.delta     = 1;
defaultOpts.tol_F     = 1.e-6;
defaultOpts.tol_opt   = 1.e-6;
defaultOpts.jacobian  = 'on'; 
defaultOpts.output    = 1;
%
if nargin == 0
   sol = defaultOpts;
   return
end 
if nargin >= 6 && ~isempty(options)
   validFields = fieldnames( defaultOpts );
   usersFields = fieldnames( options     );
   % Check if all fields in the user-provided struc are valid.
   for i = 1 : length(usersFields)
       iField = usersFields{i};
       if ~ismember( iField, validFields )
          error('%s.%s is not a valid options field.', ...
                inputname(8),iField);
       end
   end
   for i = 1 : length(validFields)
       iField = validFields{i};
       if isfield( options, iField )
          defaultOpts.(iField) = options.(iField);
       end
   end
end
% Check user-defined options 
if ~isa(defaultOpts.max_itns,'numeric')
   error('Incorrect  option.max_itns');
end
if ~isa(defaultOpts.max_feval,'numeric')
   error('Incorrect option.max_feval');
end
if ~isa(defaultOpts.delta,'numeric')
   error('Incorrect option.delta');
end
if ~isa(defaultOpts.tol_F,'numeric')
   error('Incorrect option.tol_F');
end
if ~isa(defaultOpts.tol_opt,'numeric')
   error('Incorrect option.tol_opt');
end
if ~isa(defaultOpts.jacobian,'char')
   error('Incorrect option.jacobian');
end
% 
maxit = defaultOpts.max_itns; maxnF = defaultOpts.max_feval; Delta = defaultOpts.delta;
tolf = defaultOpts.tol_F; tolg = defaultOpts.tol_opt; 
diff_jac =~ strcmp(deblank(defaultOpts.jacobian),'on');
outflag = defaultOpts.output;
%
%   Initialization
%
n = length(x); me_orig = e_i(1); mi = e_i(2); 
ind_fix = find(l==u); u_fix = u(ind_fix);              % find the fixed variables
mfix = length(ind_fix); me = me_orig+mfix; m = me+mi;  % dimensions of the problem
args = {};
if nargin>=7
   args = varargin;
end
%   Internal parameters 
infinity = 1e20;
itc = 0;  nvF = 0;  %iteration-counter, F-evaluations counter  
beta1 = 0.1; beta2 = 0.25; beta3 = 0.75; delta1 = 0.25; delta2 = 2;  %trust-region parameters
epsilon = 1e2*eps;   Deltamin1 = eps;  Deltamin2 = sqrt(eps);        %trust-region parameters
%
% Unset the bounds of the fixed variables.
%
l(ind_fix) = -infinity; u(ind_fix) = infinity;
%
if outflag >= 0
   disp(sprintf('\n')) 
   disp('******** Problem data *****************************')
%   disp(sprintf('%s %s','Problem    = ', fun))
   disp(sprintf('%s %6.0f', 'variable dimension =       ',n))
   disp(sprintf('%s %6.0f', 'number of constraints =    ',m))
   disp(sprintf('%s %6.0f', 'number of equalities =     ',me))
   disp(sprintf('%s %6.0f', 'number of inequalities =   ',mi)) 
   disp(sprintf('%s %6.0f', 'number of fixed variables =',mfix))
end
% Check on the bounds and initial guess
for i = 1 : n 
    if l(i) > u(i)
       error('ATTENTION: incorrect lower and upper bounds') 
    end
    if (x(i) < l(i) || x(i) > u(i)) 
       error('ATTENTION: the initial guess is not feasible') 
    end   
end 
% Set infinity value on the bounds;
l = max(l, -infinity); u = min(u, infinity); 
%
if diff_jac
   Cx = feval(fun,x,args{:}); 
else
   [Cx,jac] = feval(fun,x,args{:}); 
end
Fx = Cx; nvF = nvF+1; 
%
if (mi > 0 || mfix > 0)
   Fx = refF_fix_max2(Fx,x,me_orig,mi,mfix,u_fix,ind_fix);
end
Fnrm = norm(Fx); Fnrm2 = Fnrm^2; FnrmI = norm(Fx,inf);  
if mi > 0
   viol = norm(max(zeros(mi,1),Cx(me_orig+1:me_orig+mi)),inf); 
end
%
%   Jacobian evaluation (jac)
% 
if diff_jac
   if size(args) ~= 0
      jac = diffjac(x,fun,Fx,l,u,mi,me_orig,mfix,u_fix,ind_fix,varargin{:});
   else
      jac = diffjac(x,fun,Fx,l,u,mi,me_orig,mfix,u_fix,ind_fix);
   end
else 
   if (mi > 0 || mfix > 0)
      jac = refjac_fix_max2(Cx,jac,x,me_orig,mi,mfix,ind_fix);
   end
end
%   Reduced model
if mi > 0 
   [Fx,jac] = fj_rid(Fx,jac,mi); mr = size(jac); mr = mr(1); 
else 
   mr = m;
end
grad = (Fx'*jac)'; 
% 
%   Calculation of the scaling matrices  sqrt(D) and D (Dsqrt, D) 
%   and of the scaled gradient Dgrad (D.*grad)
% 
[Dsqrt,D] = scaling_matrix(x,grad,l,u);
Dsqrtgrad = Dsqrt.*grad; nDsqrtgrad = norm(Dsqrtgrad);
Dgrad = D.*grad; nDgrad = norm(Dgrad);                                     
nPgrad = norm(kk_proj(x-grad,l,u)-x); ming = min(nDgrad,nPgrad);
% 
%   Initial output
%
if outflag > 0 
   disp(sprintf('\n')) 
   disp('******** Iteration history ')
   disp(sprintf('%s %s %s %s %s %s', blanks(4),'it',blanks(3), ...
       '||F||_2', blanks(5), 'first-order', blanks(2), 'trust-region', ...
       blanks(1), 'trust-region', blanks(2), 'norm of',  blanks(6),'step' ,...
       blanks(6),'t')) 
       disp(sprintf('%s %s %s %s %s %s', blanks(25), 'optimality', blanks(5), ... 
       'radius', blanks(6), 'solution', blanks(6), 'step', blanks(5),'direction')) 
       disp(sprintf('%6.0f %s %10.4e %s %10.4e %s %4.2e',itc,blanks(2),Fnrm, ...
       blanks(4), min(nDgrad, nPgrad), blanks(4), Delta))
end
%
warning('off', 'MATLAB:nearlySingularMatrix'); 
warning('off', 'MATLAB:singularMatrix');
%   Iteration 
% 
Fnrm0 = Fnrm;
while (FnrmI > tolf && ming > tolg*sqrt(n)  && itc < maxit && nvF <= maxnF)  %WHILE#1  
% 
%  Unconstrained minimizer of min||jac*p+Fx||_2 (Newton step-pnw)
%
      if n == mr
         lastwarn('') 
         pnw = jac\(-Fx);        % square system, gaussian elimination 
         [lastmsg, lastid] = lastwarn;
         if strcmp(lastid,'MATLAB:nearlySingularMatrix') || ...
            strcmp(lastid,'MATLAB:singularMatrix')
            pnw = cod_ls(jac,-Fx);  % minimum norm solution
         end
      else
         pnw = cod_ls(jac,-Fx);  % minimum norm solution
      end
% 
%     Computation of the minimizer (pc) of the quadratic model the scaled gradient.  
% 
      vert = nDsqrtgrad^2/norm(jac*Dgrad)^2; 
      pc = -vert*Dgrad; pcv=pc;   
% 
%     Trust-region strategy 
% 
      rhof = 0; nridu = -1; alp = zeros(1,n);
      while (rhof < beta2 && Delta > Deltamin1)     % WHILE#2
            nridu = nridu+1;            
%           Computation of the generalized Cauchy point (pcg) 
            if (norm(pc)>Delta) 
               pcv = -Delta*Dgrad/nDgrad;            
            end 
            for i = 1 : n 
                if pcv(i) ~= 0  
                   alp(i) = max((l(i)-x(i))/pcv(i),(u(i)-x(i))/pcv(i)); 
                else 
                   alp(i) = infinity; 
                end 
            end         
            alpha = min(alp);
            if alpha <= 1
               pcg= alpha*pcv; 
            else
               pcg = pcv;
            end 
%        
            if norm(pnw) <= Delta
%              The Newton step is the solution of the trust-region subproblem. 
               p_tr = pnw; soltr = 'N';
            else     
%              Dogleg method 
               ngrad = norm(grad); ngrad2 = ngrad^2;
               vert = ngrad2/norm(jac*grad)^2; 
               pc_tr = -vert*grad;
               if (norm(pc_tr) >= Delta) 
%                  The Cauchy step is the solution of the trust-region subproblem 
                   p_tr = -Delta*grad/ngrad; soltr = 'C';  
               else  
                   p_tr = dogleg(pnw,pc_tr,Delta); soltr = 'D';  
               end 
            end 
%           Computation of the projected trust region step (pp_tr)
            if all(x+p_tr == kk_proj(x+p_tr,l,u))
               pp_tr = p_tr; indpr = 0;
            else
               pp_tr = kk_proj(x+p_tr,l,u)-x; indpr = 1; 
            end       
% 
%           Force condition rhoc
%  
            p1 = jac*pcg; p2 = jac*pp_tr;
            rhoc_den = grad'*pcg+0.5*norm(jac*pcg)^2;  
            rhoc = (grad'*pp_tr+0.5*norm(jac*pp_tr)^2)/rhoc_den; 
  
            if rhoc >= beta1
               ts = 0; 
            else  
               a = -0.5*norm(p1-p2)^2;
               b = -(Fx+p2)'*(p1-p2);    
               c = -grad'*pp_tr - 0.5*p2'*p2+ beta1*rhoc_den;  
               delta = b^2-4*a*c; ts = (-b+sqrt(delta))/(2*a);  
            end 
            if ts == 0  
               if indpr == 1 
                  step = 'PTR';    % projected trust-region step
               else 
                  step = ' TR';    % trust-region step
               end 
            else 
               step = ' CC';       % convex combination of TR and PTR
            end   
            p = ts*pcg+(1-ts)*pp_tr;
%
%           Accuracy requirements 
% 
            xpp = x+p;
            if diff_jac
               Cxpp = feval(fun,xpp,args{:}); 
            else
               [Cxpp,jacxpp] = feval(fun,xpp,args{:}); 
            end
            Fxpp = Cxpp; nvF = nvF+1; 
            if (mi > 0 || mfix > 0)
               Fxpp = refF_fix_max2(Fxpp,xpp,me_orig,mi,mfix,u_fix,ind_fix);
            end
            Fnrmxpp = norm(Fxpp); 
            rhof = (Fnrmxpp^2-Fnrm2)*0.5/(grad'*p+0.5*norm(jac*p)^2);  
            Deltas = Delta; Delta = min(delta1*Delta, 0.5*norm(p));   
      end   % ENDWHILE#2
% 
      if (Delta <= Deltamin1 && rhof < beta2)
         break; 
      end  
      Delta = Deltas; 
%   
%     Updating of the iterate. 
% 
      x = xpp; Cx = Cxpp; Fx = Fxpp; Fnrm = Fnrmxpp; Fnrm2 = Fnrm^2; FnrmI = norm(Fx,inf);
      if mi > 0
         viol = norm(max(zeros(mi,1),Cx(me_orig+1:me_orig+mi)),inf); 
      end
      itc = itc+1;
      if ~diff_jac
         jac = jacxpp;
      end
%
%     Jacobian evaluation (jac)
% 
      if diff_jac
         if size(args) ~= 0
            jac = diffjac(x,fun,Fx,l,u,mi,me_orig,mfix,u_fix,ind_fix,varargin{:});
         else
            jac = diffjac(x,fun,Fx,l,u,mi,me_orig,mfix,u_fix,ind_fix); 
         end
      else 
         if (mi > 0 || mfix > 0)
            jac = refjac_fix_max2(Cx,jac,x,me_orig,mi,mfix,ind_fix);
         end
      end 
%     Reduced model
      if mi > 0 
         [Fx,jac] = fj_rid(Fx,jac,mi); mr = size(jac); mr = mr(1); 
      end
      grad = (Fx'*jac)';    
% 
%     Calculation of the scaling matrices  sqrt(D) and D (Dsqrt, D) 
%     and of the scaled gradient Dgrad (Dgrad)
% 
      [Dsqrt,D] = scaling_matrix(x,grad,l,u);
      Dsqrtgrad = Dsqrt.*grad; nDsqrtgrad = norm(Dsqrtgrad);
      Dgrad = D.*grad; nDgrad = norm(Dgrad);  
      nPgrad = norm(kk_proj(x-grad,l,u)-x);  ming = min(nDgrad,nPgrad);  
% 
%     Storing and printing the iteration's summary. 
%  
      n_step=norm(p);
      if outflag > 0   
         disp(sprintf('%6.0f %s %10.4e %s %10.4e    %s %4.2e  %s %s %s %4.2e %s  %s %s %4.2e ', ...
          itc, blanks(2),Fnrm,   blanks(4), ming,  blanks(1), Delta, ...
          blanks(6), soltr, blanks(8), n_step,  blanks(4),  step, blanks(4), ts)) 
      end
      if (abs(Fnrm-Fnrm0) <= epsilon*Fnrm && Fnrm > tolf)  
	 break;
      end 
% 
%     Updating of the trust-region size 
%   
      if rhof >= beta3     
         Delta = max(Delta, max(delta2*n_step, Deltamin2)); 
      else
         Delta = max(Delta, Deltamin2); 
      end  
      Fnrm0 = Fnrm;
end       %  ENDWHILE#1
%
%  Final output 
%
sol = x;
if itc == 0
   stop = [FnrmI<tolf  ming<tolg*sqrt(n) itc >= maxit nvF>maxnF ...
 (Delta <= Deltamin1 && rhof < beta2) (abs(Fnrm-Fnrm0) <= epsilon*Fnrm && FnrmI > tolf)];
else
   stop = [FnrmI<tolf  ming<tolg*sqrt(n) itc >= maxit nvF>maxnF ...
       0 (abs(Fnrm-Fnrm0) <= epsilon*Fnrm && FnrmI > tolf)];
end

fstop = find(stop);
if fstop(1) == 1
   ierr = 0;
elseif fstop(1) == 2
       ierr = 1;
elseif fstop(1) == 3 
       ierr = 2; 
elseif fstop(1) == 4 
       ierr = 3; 
elseif fstop(1) == 5 
       ierr = 4; 
else        
    ierr = 5; 
end 
output.itns = itc;
output.feval = nvF;
output.Fnorm = Fnrm;
output.optimality = ming; 
if mi > 0
   output.violations = viol;
end  
if outflag >= 0
   disp(sprintf('\n'))
   disp('******** Final output')
   if ierr == 0
      disp(sprintf('%s %1.0f','Successful Termination. Nonlinear Residual Condition Satisfied.  Output flag= ', ierr)) 
   elseif ierr == 1
          disp(sprintf('%s %1.0f','Successful Termination. First-Order Optimality Condition Satisfied. Output Flag= ', ierr))
   elseif ierr == 2
          disp(sprintf('%s %1.0f','Maximum Number of Iterations Exceeded. Output Flag= ', ierr))   
   elseif ierr == 3
          disp(sprintf('%s %1.0f','Maximum Number of Function Evaluations Exceeded. Output Flag= ', ierr)) 
   elseif ierr == 4
          disp(sprintf('%s %1.0f','Trust-Region Radius is Less than Floating-Point Relative Accuracy. Output Flag= ', ierr))
   elseif ierr == 5
          disp(sprintf('%s %1.0f','No Improvement on Nonlinear Residual. Output Flag= ', ierr))
   end
   disp(sprintf(' '))
   disp(sprintf('%s %16.0f','Number of iterations performed: ', itc))  
   disp(sprintf('%s %2.0f','Number of function evaluations (no Jacobian): ', nvF))   
   disp(sprintf(' '))
   disp(sprintf('%s %10.5e','2-Norm of the residual F:   ', Fnrm)) 
   disp(sprintf('%s %16.5e','First-order optimality:', ming)) 
   if mi > 0
      disp(sprintf('%s %16.5e','Constraint violations: ', viol)) 
   end 
   disp('*************************************************')
end
return
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function [jac] = diffjac(x,fun,F0,l,u,mi,me,mfix,u_fix,ind_fix,varargin)
%   Compute the Jacobian matrix of function F at x by forward finite differences 
% 
%      function[jac]=diffjac(x,fun,F0,l,u,mi,me,mfix,u_fix,ind_fix,varargin); 
%  
%  inputs: 
%  x, C    = point and function. 
%  F0      = F(x), preevaluated. 
%  l,u     = constraints.
%  mi      = number of inequalities of the problem.
%  me      = number of equalities of the problem.
%  mfix    = number of fixed variables.
%  u_fix   = components of the upper bound u corresponding to the fixed variables. 
%  ind_fix = indeces corresponding to the fixed variables.  
%  varargin= additional parameters, if any, for computing C.
%
%  output  
%     jac  = approximated Jacobian matrix.  
% 
n = length(x); epsnew = sqrt(eps); jac = zeros(me+mi+mfix,n);
for j = 1 : n 
%     
%   Choice of the steplength for the forward differences 
%  
    if x(j) == 0 
       h = epsnew; 
    else 
       h = epsnew*sign(x(j))*max(abs(x(j)),norm(x,1)/n); 
    end 
    xhj = x(j)+h;  
% 
% 
    if (xhj < l(j) ||xhj > u(j)) 
% 
%      The point xhj is not feasible; backward differences is used 
% 
       h = -h; 
       xhj = x(j)+h; 
       if (xhj < l(j) || xhj > u(j)) 
          error('Function diffjac: loosing of feasibility using both backward and forward differences') 
       end 
    end 
    xh = x;  xh(j)=xhj; 
    F1=feval(fun,xh,varargin{:});  
    if mfix > 0                                      % fixed variables are present  
       F1(me+mfix+1:me+mfix+mi,1) = F1(me+1:me+mi,1); 
       for i = 1 : mfix
           F1(me+i,1) = xh(ind_fix(i))-u_fix(i);
       end
    end
    if mi > 0                                        % inequalities are present  
      F1(me+mfix+1:me+mfix+mi,1) = 0.5*max(0,F1(me+mfix+1:me+mfix+mi,1)).^2;
    end
    jac(:,j) = (F1-F0)/h;   
end 
return
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fr,jacr] = fj_rid(F0,jac,mi)
%
%  Form the reduced model when mi>0 ( a  null row in jac corresponds
%  to a zero component in F0 ).
%
%      function [Fr,jacr]=fj_rid(F0,jac,mi)
%
%  input  
%  F0   = vector function
%  jac  = jacobian matrix. 
%  mi   = number of inequalites  of the problem. 
%
%  output  
%  Fr    = reduced vector function
%  jacr  = reduced jacobian matrix. 
% 
m = size(jac,1);
z = find(F0(m-mi+1:m))+m-mi; z = z';
icol = [1:m-mi,z];
Fr = F0(icol); jacr = jac(icol,:);
return
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [Dsqrt,D] = scaling_matrix(x,grad,l,u)
%   
%   Compute the diagonal scaling matrices  (D)^{1/2}, D
% 
%      function[Dsqrt,D]=scaling_matrix(x,grad,l,u); 
% 
% input  
% x     = point. 
% grad  = gradient of f, preevaluated. 
% l,u   = constraints. 
%  
% output  
% Dsqrt =  array containing the entries of sqrt(D(x)). 
% D     =  array containing the entries of D. 
% 
infinity = 1e20;
n = length(x); Dsqrt = ones(n,1); D = ones(n,1);  
for i = 1 : n 
    if grad(i) < 0  
       if u(i) ~= infinity
          D(i) = u(i)-x(i);
          Dsqrt(i) = sqrt(D(i)); 
        end 
    else 
        if l(i) ~= -infinity 
           D(i)=x(i)-l(i); 
           Dsqrt(i) = sqrt(D(i)); 
        end 
    end 
end 
return 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function sol = dogleg(s,p,Delta) 
% 
% Compute the dogleg step. 
% 
%      function[sol]=dogleg(s,p,Delta) 
%
%  input: 
%  s     = Newton step. 
%  p     = Cauchy point. 
%  Delta = trust-region radius.
% 
%  output: 
%  sol   = dogleg step. 
% 
pnorm2 = norm(p)^2; 
a = norm(s-p)^2; b = (p'*s-pnorm2); c = pnorm2-Delta^2; 
if b == 0
    l1 = sqrt(-c/a); l2 = -l1;
else
    l1 = (-b+sign(-b)*sqrt(b^2-a*c))/a; l2 = c/(l1*a); 
end
lambda = max(l1,l2); 
sol = p+lambda*(s-p); 
return 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function px = kk_proj(x,l,u)
%  
% Projection of the point x onto the box [l,u]
%
% function px = kk_proj(x,l,u)
%
px = min(u,x); 
px = max(l,px);
return
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = cod_ls(A,b)
%
% Compute the minimum norm solution to the linear least-squares problem
% min ||Ax-b||_2by the complete orthogonal decomposition.
%              function    x=cod_ls(A,b)
% 
[U, R, V] = cod(A); %Rank decisions are made using TOL that approximates
                    %  MAX(SIZE(A))*NORM(A)*EPS.                   
[m,n] = size(A); r = length(R(1,:));
c = (b'*U)'; c = c(1:r);
y = R\c; y = [y ; zeros(n-r,1)];
x = (y'*V)';
return
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U, R, V] = cod(A, tol)
%
% Complete orthogonal decomposition.
%  
%    function  [U, R, V] = COD(A, TOL) 
%
% The decomposition has the form A = U*T*V where U and V are unitary, 
% T = [R 0; 0 0] has the same dimensions as A, and R is upper triangular 
% and nonsingular of dimension rank(A).
% Rank decisions are made using TOL, which defaults to approximately
% MAX(SIZE(A))*NORM(A)*EPS.
% By itself, COD(A, TOL) returns R.
% The function is based on the function COD contained in 
% The Matrix Computation Toolbox
% http://www.maths.manchester.ac.uk/~higham/mctoolbox/
% Slight modifications were made to fasten the original function.
% Reference:
% G.H. Golub and C.F. Van Loan, Matrix Computations, Second
% Edition, Johns Hopkins University Press, Baltimore, Maryland,
% 1989, sec. 5.4.2.
%
% input 
% A    = m-by-n matrix.
% TOL  = tolerance used to compute rank(A).
%
% output
% U, V = unitary matrices.
% R    = upper triangular and nonsingular matrix of dimension rank(A).
%
[m, n] = size(A);
% QR decomposition.
[U, R, P] = qr(A,0);    % AP = UR
%
if nargin == 1, tol = max(m,n)*eps*abs(R(1,1)); end  % |R(1,1)| approx NORM(A).
%
% Determine r = effective rank.
r = sum(abs(diag(R)) > tol);
r = r(1);             % Fix for case where R is vector.
% R = R(1:r,:);       % Throw away negligible rows (incl. all zero rows, m>n).
if r ~= n
%
%  Reduce nxr R' =  r  [L]  to lower triangular form: QR' = [Lbar].
%                 n-r  [M]                                  [0]   
   [Q, R] = trap2tri(R(1:r,:)');
   PP(P) = 1 : n;
   V = Q(:,PP);
   R = R';
else
%  V=P'; A = URV;
   I = eye(n);
   V = I(P,:);
   R = R(1:r,:);
end
if nargout <= 1, U = R; end
return
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q, T] = trap2tri(L)
%
% Unitary reduction of trapezoidal matrix to triangular form.
%
%       function [Q, T] = TRAP2TRI(L)
%
% The function is based on the function trap2tri contained in 
% The Matrix Computation Toolbox
% http://www.maths.manchester.ac.uk/~higham/mctoolbox/
% Slight modifications were made to fasten the original function.
% Reference:
% G.H. Golub and C.F. Van Loan, Matrix Computations, second edition,
% Johns Hopkins University Press, Baltimore, Maryland, 1989.
% P5.2.5, p. 220.
%
% input
% L = m-by-n lower trapezoidal  matrix with m >= n.
%
% output
% Q = unitary matrix such that QL = [T; 0]; Q is a product of 
%     Householder transformations.
% T = n-by-n and lower triangular matrix.
%         
[n, r] = size(L);
if (r > n  || norm(L-tril(L),1))
   error('Matrix must be lower trapezoidal and m-by-n with m >= n.')
end
Q = eye(n);   % To hold product of H.T.s
if r ~= n
%  Reduce nxr L =   r  [L1]  to lower triangular form: QL = [T].
%                  n-r [L2]                                 [0]
   for j = r : -1 : 1
       irow = j : n;
%  x is the vector to be reduced, which we overwrite with the H.T. vector.
       x = L(irow,j);
       x(2:r-j+1) = zeros(r-j,1);  % These elts of column left unchanged.
       s = norm(x)*(sign(x(1)) + (x(1)==0));    % Modification for sign(1)=1.
%  Nothing to do if x is zero (or x=a*e_1, but we don't check for that).
       if s ~= 0
          x(1) = x(1) + s;
          beta = s*x(1);
          xt=x'/beta;
%  Implicitly apply H.T. to pivot column.
%  L(r+1:n,j) = zeros(n-r,1); % We throw these elts away at the end.
          L(j,j) = -s;         
%  Apply H.T. to rest of matrix.
          if j > 1
             y = xt*L(irow, 1:j-1); 
             L(irow, 1:j-1) = L(irow, 1:j-1) - x*y;
          end
%  Update H.T. product. 
          y = xt*Q(irow,:); 
          for jj = j : n
              Q(jj,:) = Q(jj,:)-x(jj-j+1)*y;
          end
      end
   end
end
T = L(1:r,:);   % Rows r+1:n have been zeroed out.
return
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fx = refF_fix_max2(Fx,x,me,mi,mfix,u_fix,ind_fix)
%
% Form the residual function F when fixed variables or inequalities are present 
%  
%    function Fx=refF_fix_max2(Fx,x,me,mi,mfix,u_fix,ind_fix)     
%
% Fixed variables are transformed into equalities
% The inequalities C_I<=0 are transformed into equalities
% using the fuction t_+=0.5*max(t,0)^2
%
% input  
% Fx      = current function value.
% x       = point.
% me      = number of equalities. 
% mi      = number of inequalies.
% u_fix   = components of the upper bound u corresponding to the fixed variables. 
% ind_fix = indeces corresponding to the fixed variables.
%  
% output  
% Fx      = transformed function value.
% 
if mfix > 0                                  % fixed variables are present  
Fx(me+mfix+1:me+mfix+mi,1) = Fx(me+1:me+mi,1); 
for i = 1 : mfix
      Fx(me+i,1) = x(ind_fix(i))-u_fix(i);
end
end
if mi > 0                                    % inequalities are present  
 Fx(me+mfix+1:me+mfix+mi,1) = 0.5*max(0,Fx(me+mfix+1:me+mfix+mi,1)).^2;
end
return
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function jac = refjac_fix_max2(Cx,jac,x,me,mi,mfix,ind_fix)
%
% Form the Jacobian matrix of F when fixed variables or inequalities are present 
%
%   function jac=refjac_fix_max2(Cx,jac,x,me,mi,mfix,ind_fix)
%   
% The transpose of the gradient of the equalities 
% corresponding to the fixed variables and the transformed 
% inequalities C_I<=0 are enclosed into the analitycal Jacobian.
%
% input  
% Cx      = current function value.
% jac     = current jacobian.
% x       = point.
% me      = number of equalities.
% mi      = number of inequalies.
% ind_fix = indeces corresponding to the fixed variables.
%  
% output  
% jac     = updated jacobian
%
n = length(x);
if mfix > 0                         % fixed variables are present  
jac(me+mfix+1:me+mfix+mi,:) = jac(me+1:me+mi,:); 
jac(me+1:me+mfix,:) = zeros(mfix,n);  
for i = 1 : mfix
    jac(me+i,ind_fix(i)) = 1; 
end
end
if mi > 0                           % inequalities are present
   for i = me+mfix+1 : me+mfix+mi
       jac(i,:) = max(Cx(i-mfix),0)*jac(i,:);
   end
end
return
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

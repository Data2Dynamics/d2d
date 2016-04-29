% Wrapper for Tresnei function calling the TRESNEI.m in  with all relevant input arguments
% (C) Franz-Georg Wieland 2016
% Contact: franz-georg.wieland at mars.uni-freiburg.de 

% Necessary changes in TRESNEI.m for implementation: 
% comment out line 214: "sprintf('%s %s','Problem    = ', fun)" because sprintf 
% can't deal with function handles


function [p,exitflag,output,lambda,jac] = arTRESNEI(fun,x,l,u,arTRESoptions)

% TRESNEI allows for multidimensional constraint functions instead of
% simple constant constraints, however no support for this kind of
% constraint is implemented in D2D
% (e_i) describes the dimensions of the non-constant constraint functions
% --> set to 0 here
e_i = [0,0];

% Read out ar.optim.optimizer options and fit them to TRESNEI option structure
if isempty( arTRESoptions.MaxIter )
else options.max_itns = arTRESoptions.MaxIter;
end

if isempty( arTRESoptions.MaxFunEvals )
else options.max_feval = arTRESoptions.MaxFunEvals;
end

if isempty( arTRESoptions.TolFun )
else options.tol_F = arTRESoptions.TolFun;    
end

if isempty( arTRESoptions.Jacobian )
else options.jacobian = arTRESoptions.Jacobian;
end

% Prohibit termination of algorithm due to first order optimality criterion
% by setting its tolerance to 0
options.tol_opt = 0;


%Choose output default (documentation from TRESNEI.m)
%                         .output=0   a final summary output is printed: 
%                         .output<0   no output is displayed.
%                         .output>0   one line of summary output for each iteration is 
%                                     printed: 
options.output = 0;


% Call TRESNEI.m with relevant data, 
% x l, u are transposed because TRESNEI.m uses column notation, opposed to row notation in D2D
% Function handle is transposed for same reason
[p,ierr,output_tresnei] = TRESNEI(x',e_i,@(pTrial)transposehandle(fun, pTrial),l',u',options);


% not used
lambda = [];
jac = [];

% redefine output to fit excpected naming
output.iterations = output_tresnei.itns;                % iterations count
output.funcCount = output_tresnei.feval;                % F-evaluations count (without F-evals due to Jacobian approx
output.firstorderopt = output_tresnei.optimality;       % first-oder optimality
output.Fnorm = output_tresnei.Fnorm;                    % Fnorm (2-Norm of Function at solution)

% Translate exitflag into output Message (code taken from TRESNEI.m)
if ierr == 0
      output.ExitMessage = 'Successful Termination. Nonlinear Residual Condition Satisfied.'; 
      exitflag = 1;
   elseif ierr == 1
      output.ExitMessage ='Successful Termination. First-Order Optimality Condition Satisfied.';
      exitflag = 1;
   elseif ierr == 2
      output.ExitMessage ='Maximum Number of Iterations Exceeded.';
      exitflag = 50;
   elseif ierr == 3
      output.ExitMessage ='Maximum Number of Function Evaluations Exceeded.';
      exitflag = 51;
   elseif ierr == 4
      output.ExitMessage ='Trust-Region Radius is Less than Floating-Point Relative Accuracy.';
      exitflag = 60;
   elseif ierr == 5
      output.ExitMessage ='No Improvement on Nonlinear Residual.';
      exitflag = 1;
end
end


% Transposation of function handle because TRESNEI needs res as 1-column
% vector
function [res,sres] = transposehandle(handle, pTrial)
    
[res2,sres] = handle(pTrial);

% Transpose first output argument (not jacobian)
res = res2';
end


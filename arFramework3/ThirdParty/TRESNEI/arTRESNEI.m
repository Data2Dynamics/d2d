% Wrapper for Tresnei function calling the TRESNEI.m in  with all relevant input arguments
% 
% Necessary changes in TRESNEI.m;
% 
% comment out: "sprintf('%s %s','Problem    = ', fun)" because sprintf 
% can't deal with function handles


function [p,exitflag,output,lambda,jac] = arTRESNEI(fun,x,l,u,arTRESoptions)

% No complex constraints used in d2d thus dimensons of complex constraint
% functions (e_i) are 0
e_i = [0,0];

% Read out ar.optim.optimizer options to fit TRESNEI option structure
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

%Choose output default (documentation from TRESNEI.m)
%                         .output=0   a final summary output is printed: 
%                         .output<0   no output is displayed.
%                         .output>0   one line of summary output for each iteration is 
%                                     printed: 

options.output = 0;


% Call TRESNEI.m with relevant data
[p,ierr,output_tresnei] = TRESNEI(x',e_i,@(pTrial)transposehandle(fun, pTrial),l',u',options);


% not used
lambda = [];
jac = [];

% redefine output to fit excpected naming
output.iterations = output_tresnei.itns;
output.funcCount = output_tresnei.feval;
output.firstorderopt = output_tresnei.optimality;

% set exitflag >0, because TRESNEI.m already prints its own output

exitflag = 1;



end


% Transposation of function handle because TRESNEI needs res as 1-column
% vector
function [res,sres] = transposehandle(handle, pTrial)
    
[res2,sres] = handle(pTrial);

res = res2';
end


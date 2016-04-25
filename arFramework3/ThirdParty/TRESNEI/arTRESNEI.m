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

% Transpose first output argument (not jacobian)
res = res2';
end


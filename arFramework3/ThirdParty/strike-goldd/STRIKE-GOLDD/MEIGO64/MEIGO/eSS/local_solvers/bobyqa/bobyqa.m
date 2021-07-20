function varargout = bobyqa(fun,x0,lb,ub,varargin)
% MatLab wrapper for the mexbobyqa.F routine, which ports Powell's BOBYQA
% routine to MatLab. mexbobyqa.F needs th have been compiled via mex
% before.
%
% Input:
% fun     : objective function to be minimized
% x0      : initial guess for parameters
% lb, ub  : bounds for parameters
% options : struct with options for the algorithm:
%   Rhobeg, Rhoend, MaxFunEvals : see algorithm definition
%
% Output:
% x   : best guess for parameters
% fval: objective function at the solution, generally fval=fun(x)
% exitflag:
%   1 : The function converged to a solution x
%   0 : Number of iterations exceeded options.MaxIter or number of function
%       evaluations exceeded options.MaxFunEvals.
%   -1: The algorithm was terminated inappropriately
% output : struct with meta information:
%   funcCount   : number of function evaluations
%   algorithm   : name of the algorithm
%   t_cpu       : cpu time
%
% Currently has 2nd mode with 0 input arguments returning a function handle
% used by funHandleWrap.m to compute objective function values. This is
% because via mexCallMATLAB only non-anonymous functions (e.g. defined in a
% file) or global function handles can be called, but pesto uses anonymous
% function handles. ATTENTION: This might be non-thread-safe.
%
% The current interface to Powell'S BOBYQA routine is rather rudimentary,
% so feel free to improve on this.

%% save objective funtion

% sorry for being global
global bobyqafun_sfbg;

if nargin == 0
    varargout{1} = bobyqafun_sfbg;
    return;
end

%% main part

if ~exist('mexbobyqa', 'file')
    error(sprintf(['The mexbobyqa file does not exist. Please compile it first\n',...
		'by running the file compile_mexbobyqa.m in the bobyqa directory.']));
end

bobyqafun_sfbg = fun;

% check for options
if (nargin > 4)
    options = varargin{1};
else
    options = struct();
end

% interpret parameters

N       = length(x0);

if (isfield(options,'Npt') && ~isempty(options.Npt))
    NPT		= options.Npt;
else
    NPT  = 2*N+1;
end

X       = x0;
LB      = lb;
UB      = ub;

CALFUN  = func2str(fun);

if (isfield(options,'Rhobeg') && ~isempty(options.Rhobeg))
    RHOBEG		= options.Rhobeg;
else
    RHOBEG  = 1e-1;
end

if (isfield(options,'Rhoend') && ~isempty(options.Rhoend))
    RHOEND		= options.Rhoend;
else
    RHOEND  = 1e-8;
end

if (isfield(options,'MaxFunEvals') && ~isempty(options.MaxFunEvals))
    MAXFUN		= options.MaxFunEvals;
else
    MAXFUN  = 1000*N;
end

IPRINT = 0; % no output

% track time
starttime = cputime;

% do optimization
[ x,fval,feval ] = mexbobyqa(N,NPT,X,LB,UB,CALFUN,RHOBEG,RHOEND,MAXFUN,IPRINT);

% meta information
output.funcCount = feval;
output.algorithm = 'BOBYQA';
output.t_cpu = cputime - starttime;
exitflag = 1; % extract it from fortran if you need it

% return
varargout{1} = x;
varargout{2} = fval;
varargout{3} = exitflag;
varargout{4} = output;
end
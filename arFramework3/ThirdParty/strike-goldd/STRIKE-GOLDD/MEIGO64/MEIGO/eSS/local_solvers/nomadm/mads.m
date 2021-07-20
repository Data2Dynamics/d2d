function [BestF,BestI,RunData,varargout] = mads(Problem,iterate0,varargin)
%MADS  Solver for nonlinear and mixed variable constrained optimization
%
%   Syntax:
%      [BESTF,BESTI,RUNDATA] = mads(PROBLEM,ITERATE0)
%      [BESTF,BESTI,RUNDATA] = mads(PROBLEM,ITERATE0,OPTIONS)
%      [BESTF,BESTI,RUNDATA] = mads(PROBLEM,ITERATE0,OPTIONS,RUNDATA)
%      [BESTF,BESTI,RUNDATA,CACHE] = mads(PROBLEM,ITERATE0)
%      [BESTF,BESTI,RUNDATA,CACHE] = mads(PROBLEM,ITERATE0,OPTIONS)
%      [BESTF,BESTI,RUNDATA,CACHE] = mads(PROBLEM,ITERATE0,OPTIONS,RUNDATA)
%
%   Description:
%      MADS is an optimization algorithm for constrained nonlinear and mixed
%      variable programming problems.  It employs the derivative-free class of
%      mesh-adaptive direct search (MADS) filter algorithms, which are a
%      generalization of the class of Generalized Pattern Search (GPS) methods.
%      Derivatives are not necessary to run the algorithm, but if available can
%      be used to make the algorithm faster and more efficient in finding a
%      solution.
%
%      All the input and output variables are structures, each of which
%      contains many parameters.  The variables BESTF and BESTI, contain the
%      data describing the best feasible point and least infeasible point found
%      by the algorithm.  RUNDATA contains statistics on the MADS run, such as
%      final mesh size, number of iterations, number of function and gradient
%      evaluations, and CPU time used.  CACHE contains data associated with
%      every point evaluated during the run.
%
%      The input variable PROBLEM contains variables that describe the
%      optimization problem to be solved.  In particular, it contains the names
%      of all the files associated with the problem, such as the functions
%      file, the linear constraints file, the discrete neighbors file for MVP
%      problems, and the initial points file.  ITERATE0 is a structure
%      describing all the initial iterates.  OPTIONS contains all user options
%      and settings to be used during the run.  If not used as an input
%      argument, MADS uses default settings, as specified in the MADS_DEFAULTS
%      function.  Using RUNDATA as an input variable is only done within the
%      NOMADm GUI setup, when resuming a previously stopped run.  It is not
%      recommended for running in batch mode.
%
%      The user is referred to the User's Guide for more detailed help in
%      understanding the input and output variables.
%
%   See also MADS_BATCH, MADS_DEFAULTS, NOMADM

%*******************************************************************************
%   Copyright (c) 2001-2005 by Mark A. Abramson
%
%   This file is part of the NOMADm software package.
%
%   NOMADm is free software; you can redistribute it and/or modify it under the
%   terms of the GNU General Public License as published by the Free Software
%   Foundation; either version 2 of the License, or (at your option) any later
%   version.
%
%   NOMADm is distributed in the hope that it will be useful, but WITHOUT ANY
%   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
%   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
%   details.
%
%   You should have received a copy of the GNU General Public License along with
%   NOMADm; if not, write to the Free Software Foundation, Inc., 
%   59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% ------------------------------------------------------------------------------
%   Originally created, 2001.
%   Last modified, 2 February 2005
%
%   Author information:
%   Mark A. Abramson, LtCol, USAF, PhD
%   Air Force Institute of Technology
%   Department of Mathematics and Statistics
%   2950 Hobson Way
%   Wright-Patterson AFB, OH 45433
%   (937) 255-3636 x4524
%   Mark.Abramson@afit.edu
%*******************************************************************************

%*******************************************************************************
% mads: Runs the Mesh-Adpative Direct Search (MADS) Filter Algorithm for
%       solving constrained Mixed Variable Programming (MVP) problems.
% ------------------------------------------------------------------------------
% CALLED FUNCTIONS:
%  processInput           = Initialize and validate input data and CACHE
%    mads_defaults        =   Retrieves MADS parameter default settings
%    checkErrors          =   Error check the input data
%    createCache          =   Create and initialize the Cache
%    createFilter         =   Create and initialize a filter
%    makeFeasible         =   Convert infeasible starting point to feasible
%  processOutput          = Process output and delete temporary files/variables
%    plotHistory          =   Plot objective value vs number of function evals
%    closeWorkspace       =   Close appdata and delete temporary files
%  search                 = Perform MADS Search step
%    lhs                  =   Perform Latin Hypercube Search on the MADS mesh
%    ga                   =   Perform Genetic Algorithm
%      setPenaltyTerm     =     Set penalty parameter for constrained problem
%    recalSurrogate       =   Construct or recalibrate a surrogate
%    optimizeSurrogate    =   Solve a surrogate optimization problem
%  poll                   = Perform MADS Poll step
%    standDirections      =   Retrieve standard mesh directions
%      activeConstraints  =     Identify e-active (linear) constraints
%      getTangentCone     =     Compute generators for the tangent cone
%        removeRedundancy =       Remove redundant active constraints
%      scaleDirections    =     Scale the directions used in the Poll step
%    gradDirections       =   Retrieve gradient-pruned mesh directions
%      dVector            =     Compute the appropriate descent direction
%    madsDirections       =   Retrieve MADS mesh directions
%    getPollOrder         =   Set the order in which iterates are polled
%  mvpPoll                = Perform discrete neighbor and extended polling
%    < User N >           =   User-defined problem discrete neighbor file
%  update                 = Update MADS parameters and run statistics
%    getPollCenter        =   Determine the next poll center
%  updateOmega            = Update Omega and scale factors
%    getScaleFactors      =   Compute scale factors, based on initial data
%    < User O >           =   User-defined problem Omega file
%  evalPointSet           = Evaluate functions at multiple trial points
%    inOmega              =   Test if feasible w.r.t. linear constraints
%    isCacheHit           =   Test if point is already in the Cache
%    evalFunction         =   Evaluate problem functions at a given point
%      < User F >         =     User-defined problem Function file
%    updateCache          =   Update the Cache with the current iterate
%    updateFilter         =   Update the filter and solutions vectors
%      dominates          =     Tests if one iterate dominates another
%      plotFilter         =     Update the current real-time filter plot
%  evalRSPointSet         = Evalaute stochastic functions at multiple points
%    RS                   =   Perform ranking and selection procedure
%  terminate              = Tests if termination criteria are satisfied
% ------------------------------------------------------------------------------
% USER FUNCTION FILES needed for Optimization Problem: example
%   example.m        = returns f(x,p), c(x,p), f'(x,p), c'(x,p)
%   example_Omega.m  = returns A,l,u,plist, given parameter p
%   example_x0.m     = returns initial iterate(s)
%   example_N.m      = returns the set N of neighbor points 
%   example_Param.m  = returns a structure Param of user parameters
%   <Search file>    = optional function specifying Search procedure
% ------------------------------------------------------------------------------
% VARIABLES (only for the mads function):
%  BestF        = final best feasible solution found
%  BestI        = final least infeasible solution found
%  RunData      = statistics measuring algorithm performance
%    .delta     =   mesh size
%    .nIter     =   cumulative number of iterations
%    .nFunc     =   cumulative number of function evaluations
%    .time      =   cumulative amount of CPU time expired
%    .nFails    =   number of consecutive Poll failures
%    .pcenter   =   current poll center
%    .stopRun   =   flag for stopping run at the next iteration
%  Cache        = structure containing previously evaluated iterates
%    .Filter    =   structure containing filter data
%  Problem      = structure containing optimization problem data
%    .nameCache =   name of the base workspace Cache variable
%    .isMVP     =   flag indicating a mixed variable problem
%  iterate0     = initial iterate
%  Options      = structure containing SEARCH and POLL parameters
%    .Term      =   criteria for terminating the MADS function
%  terminate    = inline function used to test stopping criteria
%  success      = current iteration success or failure
%*******************************************************************************

% Set up and initialize inline termination function
terminate = inline(['R.delta  <= T.delta  || R.nIter >= T.nIter || ', ...
        'R.nFunc  >= T.nFunc  || R.time  >= T.time  || ', ...
        'R.nFails >= T.nFails || R.stopRun'], 'T', 'R');

% Initialize and verify validity of program variables
[Problem,Options,RunData] = processInput(Problem,iterate0,varargin{:});
tic;

try
    
    % Main loop: Search, Poll, and Update parameters
    while(~terminate(Options.Term,RunData))
        Cache = getappdata(0,Problem.nameCache);
        [success,TempFilter] = search(Problem,Options,RunData,Cache.Filter);
        Cache = getappdata(0,Problem.nameCache);
        Cache.Filter = TempFilter;
        setappdata(0,Problem.nameCache,Cache);
        if (~success)
            [success,TempFilter,RunData] = poll(Problem,Options,RunData,...
                RunData.pcenter,Cache.Filter);
            Cache = getappdata(0,Problem.nameCache);
            Cache.Filter = TempFilter;
            setappdata(0,Problem.nameCache,Cache);
        end
        if (~success && Problem.isMVP)
            success        = mvpPoll(Problem,Options,RunData);
        end
        [Problem,RunData] = update(Problem,Options,RunData,success);
    end
    [Cache,BestF,BestI] = processOutput(Problem,Options);
    varargout = {Cache};
catch
    Cache = closeWorkspace(Problem);
    varargout = {Cache};
    rethrow(lasterror);
end
return;

%*******************************************************************************
% FUNCTIONS FOR PROCESSING INPUT DATA
%*******************************************************************************
%*******************************************************************************
% processInput:  Initializes, validates, and processes input data
% ------------------------------------------------------------------------------
% Called by: mads
% Calls:     madsDefaults, checkErrors,  updateOmega,   createCache,
%            evalPointSet, makeFeasible, getPollCenter, update
% VARIABLES:
%  Problem          = structure containing optimization problem data
%    .File          =   structure of file names
%      .F           =     name of functions file
%      .O           =     name of Omega file
%      .N           =     name of discrete neighbor file
%      .C           =     name of pre-existing Cache file
%      .P           =     name of parameter file
%    .nameCache     =   name of the base workspace Cache variable
%    .fType         =   type of functions file (C=C, F=FORTRAN, M=MATLAB)
%    .Omega         =   structure defining linear constraints l <=Ax <= u
%    .iterate0      =   initial iterates
%    .isMVP         =   flag indicating an MVP problem
%    .adaptOmega    =   flag for updating an MVP's Omega parameters
%    .maxNp         =   maximum number of categorical variables
%    .maxNx         =   maximum number of continuous variables
%    .maxNc         =   maximum number of nonlinear constraints
%    .pollBasis     =   basis used to construct custom Poll directions
%  Options          = structure of MADS parameters
%    .tolBind       =   tolerance for flagging active linear constraints
%    .pollStrategy  =   string code indicating choice of Poll strategy
%    .pollOrder     =   string code indicating choice of Polling order
%    .pollCenter    =   number of the filter point that is Poll center
%    .delta0        =   initial mesh size
%    .deltaMax      =   maximum mesh size
%    .meshRefine    =   mesh refinement factor
%    .meshCoarsen   =   mesh coarsening factor
%    .hmin          =   minimum constraint violation of an infeasible point
%    .hmax          =   maximum allowable constraint violation
%    .Term          =   criteria for terminating the MADS function
%      .delta       =     lowest allowed mesh size
%      .nIter       =     maximum number of iterations
%      .nFunc       =     maximum number of function calls
%      .time        =     maximum allowable CPU time
%      .nFails      =     maximum number of consecutive Poll failures
%    .ePollTriggerF =   objective function extended poll trigger
%    .ePollTriggerH =   constraint violation function extended poll trigger
%    .nSearches     =   number of different SEARCH types to be used
%    .Search(n)     =   structure containing parameters for each SEARCH
%      .type        =     string code indicating choice of SEARCH strategy
%      .label       =     long text label of SEARCH strategy
%      .nIter       =     number of iterations to perform the SEARCH
%      .nPoints     =     maximum number of points evaluated in SEARCH step
%      .file        =     name of optional file containing the user SEARCH
%    .computeGrad   =   flag for computing any available gradient
%    .accelerate    =   flag for accelerating mesh refinement
%    .useFilter     =   turns on/off the filter for nonlinear constraints
%    .runStochastic =   flag to run as a stochastic optimization problem
%    .plotFilter    =   turns on/off real-time filter plot
%    .fplothandle   =   handle of the filter plot
%    .loadCache     =   flag for loading Cache of iterates
%  RunData          = structure of MADS run statistics and parameters
%    .porder        =   fixed Poll order
%    .scale         =   scale factors for scaling mesh directions
%    .stopRun       =   flag for stopping MADS run immediately
%  iterate0         = initial iterates
%  Cache            = structure of all previously evaluated iterates
%    .iterate       =   vector of all previously computed iterates
%    .Filter        =   structure of filter data
%      .plot        =     turns on/off filter plot
%      .plothandle  =     handle for the axes of the plot
%    .size          =   number of points in the Cache
%  Defaults         = defaults for all MADS parameters
%    .Options       =   Default values for Options
%    .Types         =   list of strings corresponding to possible choices
%      .poll        =     possible poll strategies
%      .pollOrder   =   possible poll order strategies
%  n0               = number of initial iterates
%  ind              = index that identifies which error was flagged
%  TempFilter       = used to reconstruct a filter from a previous run
%    .plothandle    = temporary storage for Cache.Filter.plothandle
%  Param            = user-defined parameters
%    .maxNp         =   maximum number of categorical variables
%    .maxNx         =   maximum number of continuous variables
%    .pollBasis     =   basis used to construct custom Poll directions
%    .pollOrder     =   fixed Poll order
%  iterate          = initial iterate after being processed
%*******************************************************************************
function [Problem,Options,RunData] = processInput(Problem,iterate0,varargin)

% Set Options, depending on the number of input arguments
Defaults = mads_defaults('Truth');
switch nargin
    case {2}
        Options = Defaults.Options;
    case {3,4}
        Options = deal(varargin{1});
        Options.tolBind = Defaults.Options.tolBind;
    otherwise
        error('mads:input','Incorrect number of input arguments (processInput).');
end

% Process initial points into the proper format
n0 = length(iterate0);
if isstruct(iterate0(1))
    [Problem.iterate0(1:n0).x] = deal(iterate0.x);
    [Problem.iterate0(1:n0).p] = deal(iterate0.p);
else
    [Problem.iterate0(1:n0).x] = deal(iterate0);
    [Problem.iterate0(1:n0).p] = deal({});
end
for k = 1:n0
    Problem.iterate0(k).n = length(iterate0(k).x);
end

% Check input data for errors
Problem.isMVP       = ~~exist(Problem.File.N,'file');
Options.computeGrad = strncmpi(Options.pollStrategy,'Gradient',8);
checkErrors(Problem,Options,Defaults.Types);

% Add a bogus empty Search to the end of the list, for convenience
if (Options.nSearches > 0)
    Options.Search(end+1) = struct('type','None','label','None', ...
        'nIter',Inf,'nPoints',0,'file','',...
        'local',0,'merit',0);
    Options.nSearches = Options.nSearches + 1;
end

% Adjust choices for specific Poll choices
if (Options.computeGrad && Problem.fType == 'M' && ...
        abs(nargout(Problem.File.F)) <= 2)
    Options.computeGrad  = 0;
    Options.pollStrategy = 'Standard_2n';
    warning('Derivatives are not available; using standard 2n directions');
end
if strncmpi(Options.pollStrategy,'MADS',4)
    Options.meshCoarsen = 4.0;
    Options.meshRefine  = 0.25;
    Options.accelerate  = 0;
end

% Change options if a filter is not used
if (~Options.useFilter)
    Options.pollCenter = 0;
    Options.plotFilter = 0;
    for k = 1:Options.nSearches
        if (Options.Search(k).type(2) == 'P')
            Options.Search(k).type = 'None';
        end
    end
end

% Clear any previous appdata
if isappdata(0,Problem.nameCache)
    rmappdata(0,Problem.nameCache);
end
if isappdata(0,'SUR')
    rmappdata(0,'SUR');
end

% Initialize Cache or load a pre-existing one
Cache = createCache(Options,Problem.File.C);
setappdata(0,Problem.nameCache,Cache);

%Initialize R&S Parameters
if Options.runStochastic
    RunData.RS           = Options.RS;
    RunData.RS.iz        = RunData.RS.iz_const;
    RunData.RS.alpha     = RunData.RS.alpha_const;
    RunData.RS.F         = [];
    RunData.RS.nFuncLeft = Options.Term.nFunc;
else
    RunData.RS = [];
end

% Set memory allocation variables and user-defined parameters
Problem.maxNp = 2*length(Problem.iterate0(1).p);
Problem.maxNx = (Problem.iterate0(1).n)*(2*Problem.isMVP + ~Problem.isMVP);
if (exist(Problem.File.P,'file') == 2)
    Param = feval(Problem.File.P);
    setappdata(0,'PARAM',Param);
    if (isfield(Param,'MaxNp')),     Problem.maxNp     = Param.maxNp; end
    if (isfield(Param,'MaxNx')),     Problem.maxNx     = Param.maxNx; end
    if (isfield(Param,'PollBasis')), Problem.pollBasis = Param.pollBasis; end
    if (isfield(Param,'PollOrder')), RunData.porder    = Param.pollOrder; end
end

% Flag Custom Poll Order or Directions error
if strncmp(Options.pollStrategy,'Custom',6) && ~isfield(Problem,'PollBasis')
    error('mads:user','User-provided Poll Basis not found (processInput).');
end
if strcmp(Options.pollOrder,'Custom') && ~isfield(RunData,'porder')
    error('mads:user','User-provided Poll Order not found (processInput).');
end

% Initialize Omega and process initial iterates
[Problem.Omega,RunData.scale] = updateOmega(Problem,Options, ...
    Problem.iterate0(1));
if exist(Problem.File.O,'file')
    Problem.adaptOmega = nargin(Problem.File.O) > 1;
else
    Problem.adaptOmega = 0;
end

% Process initial iterates
[iterate,success,TempFilter,RunData.RS] = ...
    evalPointSet(Problem,Problem.iterate0,1,Options,Cache.Filter,RunData.RS,2);
Cache = getappdata(0,Problem.nameCache);
Cache.Filter = TempFilter;
setappdata(0,Problem.nameCache,Cache);

% Create a valid Cache point if none exists
if (~success || all(~isfinite([Cache.iterate(1:Cache.size).f])))
    iterate = makeFeasible(Problem.iterate0(1),Problem.Omega);
    [iterate,success,TempFilter,RunData.RS] = ...
        evalPointSet(Problem,iterate,1,Options,Cache.Filter,RunData.RS,2);
    Cache = getappdata(0,Problem.nameCache);
    Cache.Filter = TempFilter;
    setappdata(0,Problem.nameCache,Cache);
end
Problem.maxNc = length(iterate(1).c)*(2*Problem.isMVP + ~Problem.isMVP);

% Initialize RunData parameters
if (nargin >= 4)                      % For resuming previous run
    RunData = varargin{2};
    RunData.stopRun = 0;
else
    [Problem,RunData] = update(Problem,Options,RunData,-1);
end

return;

%*******************************************************************************
% checkErrors:  Check for errors in the input data.
% ------------------------------------------------------------------------------
% Called by: processInput
% VARIABLES:
%  Problem          = structure containing optimization problem data
%    .File          =   structure of file names
%      .F           =     name of functions file
%      .O           =     name of Omega file
%      .N           =     name of discrete neighbor file
%      .C           =     name of pre-existing Cache file
%      .P           =     name of parameter file
%    .nameCache     =   name of the base workspace Cache variable
%    .fType         =   type of functions file (C=C, F=FORTRAN, M=MATLAB)
%    .Omega         =   structure defining linear constraints l <=Ax <= u
%    .iterate0      =   initial iterates
%    .isMVP         =   flag indicating an MVP problem
%    .adaptOmega    =   flag for updating an MVP's Omega parameters
%    .maxNp         =   maximum number of categorical variables
%    .maxNx         =   maximum number of continuous variables
%    .maxNc         =   maximum number of nonlinear constraints
%    .pollBasis     =   basis used to construct custom Poll directions
%  Options          = structure of MADS parameters
%    .tolBind       =   tolerance for flagging active linear constraints
%    .pollStrategy  =   string code indicating choice of Poll strategy
%    .pollOrder     =   string code indicating choice of Polling order
%    .pollCenter    =   number of the filter point that is Poll center
%    .delta0        =   initial mesh size
%    .deltaMax      =   maximum mesh size
%    .meshRefine    =   mesh refinement factor
%    .meshCoarsen   =   mesh coarsening factor
%    .hmin          =   minimum constraint violation of an infeasible point
%    .hmax          =   maximum allowable constraint violation
%    .Term          =   criteria for terminating the MADS function
%      .delta       =     lowest allowed mesh size
%      .nIter       =     maximum number of iterations
%      .nFunc       =     maximum number of function calls
%      .time        =     maximum allowable CPU time
%      .nFails      =     maximum number of consecutive Poll failures
%    .ePollTriggerF =   objective function extended poll trigger
%    .ePollTriggerH =   constraint violation function extended poll trigger
%    .nSearches     =   number of different SEARCH types to be used
%    .Search(n)     =   structure containing parameters for each SEARCH
%      .type        =     string code indicating choice of SEARCH strategy
%      .label       =     long text label of SEARCH strategy
%      .nIter       =     number of iterations to perform the SEARCH
%      .nPoints     =     maximum number of points evaluated in SEARCH step
%      .file        =     name of optional file containing the user SEARCH
%    .computeGrad   =   flag for computing any available gradient
%    .accelerate    =   flag for accelerating mesh refinement
%    .useFilter     =   turns on/off the filter for nonlinear constraints
%    .plotFilter    =   turns on/off real-time filter plot
%    .fplothandle   =   handle of the filter plot
%    .loadCache     =   flag for loading Cache of iterates 
%  indNone          = indices of Searches that are of type "None" 
%  errorMsg         = cell array of error messages
%  errorCheck       = vector of flags that track input data errors
%*******************************************************************************
function checkErrors(Problem,Options,Types)

% Check input parameters for errors and display, if any
errorMsg{1}  = 'Required Functions file not found';
errorMsg{2}  = 'Required Omega file not found';
errorMsg{3}  = 'Required discrete neighbors file not found';
errorMsg{4}  = 'MVP initial iterate is missing categorical variables';
errorMsg{5}  = 'Invalid choice of Search Method';
errorMsg{6}  = 'Invalid choice of Poll Directions';
errorMsg{7}  = 'Invalid choice of Polling Order';
errorMsg{8}  = 'Invalid choice of Poll Center';
errorMsg{9}  = 'Termination criteria must be positive numbers or infinity';
errorMsg{10} = 'Termination criteria cannot all be infinite';
errorMsg{11} = 'Initial Mesh Size must be a strictly positive number';
errorMsg{12} = 'Maximum Mesh Size must be greater than Initial Mesh Size';
errorMsg{13} = 'Mesh Refinement Factor must be a rational number in (0,1)';
errorMsg{14} = 'Mesh Coarsening Factor must be a rational number > 1';
errorMsg{15} = 'Min Filter Constraint Violation must be nonnegative';
errorMsg{16} = 'Max Filter Constraint Violation must be positive';
errorMsg{17} = 'Filter Constraint Violation Values must have Min < Max';
errorMsg{18} = 'Extended Poll Triggers must be positive numbers';
errorMsg{19} = 'Number of R&S samples must be positive';
errorMsg{20} = 'Initial Indifference Zone parameter must be > 0';
errorMsg{21} = 'Initial Alpha parameter must be a number in (0,1)';
errorMsg{22} = 'R&S Decay Factors must be numbers in (0,1)';
errorCheck(1) =~exist([Problem.File.F, '.', lower(Problem.fType)],'file');
errorCheck(2) = exist(Problem.File.N,'file') && ~exist(Problem.File.O,'file');
errorCheck(3) =~exist(Problem.File.N,'file') && ~isempty(Problem.iterate0(1).p);
errorCheck(4) = exist(Problem.File.N,'file') &&  isempty(Problem.iterate0(1).p);
errorCheck(5) = Options.nSearches > 0 && ...
    any(~ismember({Options.Search(:).type}, Types.search));
errorCheck(6)  = ~ismember(Options.pollStrategy, Types.poll);
errorCheck(7)  = ~ismember(Options.pollOrder,    Types.pollOrder);
errorCheck(8)  = ~isfinite(Options.pollCenter) || Options.pollCenter < 0;
errorCheck(9)  = Options.Term.delta <= 0       || Options.Term.nIter <= 0 ...
    || Options.Term.nFunc <= 0       || Options.Term.time  <= 0 ...
    || isnan(Options.Term.delta)     || isnan(Options.Term.nIter) ...
    || isnan(Options.Term.nFunc)     || isnan(Options.Term.time);
errorCheck(10) = isinf(Options.Term.delta)     && isinf(Options.Term.nIter) ...
    && isinf(Options.Term.nFunc)     && isinf(Options.Term.time);
errorCheck(11) = ~isfinite(Options.delta0)     || Options.delta0 <= 0;
errorCheck(12) = Options.delta0 > Options.deltaMax;
errorCheck(13) = ~isfinite(Options.meshRefine) ...
    || Options.meshRefine <= 0 || Options.meshRefine >= 1;
errorChech(14) = ~isfinite(Options.meshCoarsen) || Options.meshCoarsen < 1;
errorCheck(15) = ~isfinite(Options.hmin)        || Options.hmin <  0;
errorCheck(16) = ~isfinite(Options.hmax)        || Options.hmax <= 0;
errorCheck(17) = Options.hmax <= Options.hmin;
errorCheck(18) = isnan(Options.ePollTriggerF)  || Options.ePollTriggerF <= 0 ...
    || isnan(Options.ePollTriggerH)  || Options.ePollTriggerH <= 0;
if Options.runStochastic
    errorCheck(19) = ~isfinite(Options.RS.s0)       || Options.RS.s0 <= 0;
    errorCheck(20) = ~isfinite(Options.RS.iz_const) || Options.RS.iz_const <= 0;
    errorCheck(21) = ~isfinite(Options.RS.alpha_const) ...
        || Options.RS.alpha_const <= 0 || Options.RS.alpha_const >= 1;
    errorCheck(22) = ~isfinite(Options.RS.iz_rho) ...
        || ~isfinite(Options.RS.alpha_rho) ...
        || Options.RS.iz_rho <= 0  || Options.RS.alpha_rho <= 0 ...
        || Options.RS.iz_rho >= 1  || Options.RS.alpha_rho >= 1;
end
ind = find(errorCheck);
if ~isempty(ind),
    error('mads:input',[errorMsg{ind(1)},' (checkErrors).']);
end

return

%*******************************************************************************
% createCache:  Create and initialize the Cache.
% ------------------------------------------------------------------------------
% Called by: processInput
% Calls:     createFilter
% VARIABLES:
%  Cache            = structure containing a Cache of iterates
%    .tol           =   tolerance for an iterate being in the cache
%    .iterate       =   vector of iterates
%      .n           =     size of the continuous variable space
%      .x(nc)       =     values of continuous variables 
%      .p(nd)       =     values of categorical variables
%      .f           =     objective function value (f-value)
%      .c(m)        =     constraint function values (c(x) <= 0)
%      .h           =     constraint violation function value (h-value)
%      .gradf(nc)   =     partial derivatives of f (if available)
%      .gradc(m,nc) =     partial derivatives of c (if available)
%      .gradh       =     partial derivatives of h (if available)
%    .size          =   number of iterates
%    .xnorm         =   vector of inf-norms of the iterates
%    .pID           =   vector of unique IDs for categorical variables
%    .isSurrPoint   =   vector of flags of points used in surrogate recal
%    .bfp           =   vector of best cumulative feasible point f-values
%    .lip           =   matrix of least infeasible point f- and h-values
%    .nHits         =   number of Cache hits
%    .Filter        =   structure containing filter data
%  Options          = structure containing SEARCH and POLL parameters
%    .computeGrad   =   logical indicating gradients are to be computed
%    .hmin          =   minimum constraint violation of an infeasible point
%    .hmax          =   maximum allowable constraint violation
%    .plotFilter    =   logical for plotting real-time filter
%    .fplothandle   =   handle for the real-time filter plot
%    .tolCache      =   tolerance for a point to be a Cache point
%    .CacheFile     =   name of file containing a pre-existing Cache
%*******************************************************************************
function Cache = createCache(Options,CacheFile)

% A Cache file exists and will be used
if (Options.loadCache && exist(CacheFile,'file'))
    load(CacheFile,'Cache');
    TempFilter = createFilter(0.1,1,0,1,Options.fplothandle);
    Cache.Filter.plothandle = TempFilter.plothandle;
    if (Cache.Filter.plot)
        plotFilter(Cache.Filter,Cache.Filter.hmax,Cache.iterate,10);
    end
    
    % A new Cache is constructed
else
    Cache = [];
    Cache.tol                 = Options.tolCache;
    Cache.iterate(1).x        = [];
    Cache.iterate(1).p        = {};
    Cache.iterate(1).n        = 0;
    Cache.iterate(1).f        = [];
    Cache.iterate(1).c        = [];
    Cache.iterate(1).h        = [];
    if (Options.computeGrad)
        Cache.iterate(1).gradf = [];
        Cache.iterate(1).gradc = [];
        Cache.iterate(1).gradh = [];
    end
    Cache.size                = 0;
    Cache.xnorm               = [];
    Cache.pID                 = [];
    Cache.isSurrPoint         = [];
    Cache.nFunc               = [];
    Cache.bfp                 = [];
    Cache.lip                 = [];
    Cache.nHits               = 0;
    Cache.iterate(1)          = [];            % This is a key statement!
    Cache.Filter = createFilter(Options.hmin,Options.hmax,0,...
        Options.plotFilter,Options.fplothandle);
end
return

%*******************************************************************************
% createFilter:  Create and initialize a filter.
% ------------------------------------------------------------------------------
% Called by: processInput, mvpPoll
% VARIABLES:
%  Filter        = structure containing the filter
%    .feasible   =   indices of feasible iterates, sorted by f-value
%    .F          =   indices of iterates in the filter, sorted by h-value
%    .hmin       =   minimum allowable h-value to be in Filter.F
%    .hmax       =   maximum allowable h-value to be in Filter.F
%    .strict     =   flag indicating that hmax is strictly enforced
%    .plot       =   flag for plotting real-time filter
%    .plothandle =   handle of the real-time filter plot
%  hmin          =   input value for Filter.hmin
%  hmax          =   input value for Filter.hmax
%  strict        =   input value for Filter.strict
%  plot          =   input value for Filter.plot
%  plothandle    =   input value for Filter.plothandle
%*******************************************************************************
function Filter   = createFilter(hmin,hmax,strict,fplot,plothandle)
if (hmin >= hmax)
    error('mads:filter',['Attempt to create filter with ',...
            'incompatible hmin/hmax (createFilter).']);
end
Filter.hmin       = hmin;
Filter.hmax       = hmax;
Filter.strict     = strict;
Filter.plot       = fplot;
Filter.F          = [];
Filter.feasible   = [];
if (isempty(plothandle) && fplot)
    figure; plothandle = gca;
    set([plothandle; get(plothandle,'Children')],'Visible','off');
end
Filter.plothandle = plothandle;
return

%*******************************************************************************
% makeFeasible:  Projects an infeasible iterate (w.r.t. Omega) into Omega.
% ------------------------------------------------------------------------------
% Called by: processInput
% Calls:     inOmega
% VARIABLES:
%  iterate    = current iterate that will be made feasible
%    .x       =   continuous variables
%    .p       =   categorical variables
%  Omega      = structure defining linear constraints, l <= A*x <= u
%    .A       =   matrix of linear coefficients
%    .l       =   vector of lower bounds
%    .u       =   vector of upper bounds
%    .plist   =   list of allowable discrete variable values
%  pass       = flag indicating if categorical variables are feasible
%  [m,n]      = dimensions of linear constraint matrix Omega.A
%  ind        = indices of infeasible lower and upper bounds
%  c,A,b      = matrices that define the LP, min {c'x: Ax <= b}
%  OptOptions = structure used in LP solver to set optimization options
%  dist       = closest distance of iterate to Omega
%  err        = error flag for LP solver
%*******************************************************************************
function iterate = makeFeasible(iterate,Omega)

% Adjust categorical variables
for k = 1:length(iterate.p)
    if (ischar(iterate.p{k}))
        pass = any(strcmp(iterate.p{k},Omega.plist{k}));
    else
        pass = ~isempty(find([Omega.plist{k}{:}] == iterate.p{k},1));
    end
    if (~pass)
        iterate.p{k} = Omega.plist{k}{1};
    end
end

% Adjust bounds if necessary
[m,n] = size(Omega.A);
if n <= m && all(all(Omega.A(1:n,:) == eye(n)))
    ind = find(iterate.x - Omega.l(1:n) < 0.0);
    iterate.x(ind) = deal(Omega.l(ind));
    ind = find(iterate.x - Omega.u(1:n) > 0.0);
    iterate.x(ind) = deal(Omega.u(ind));
    if inOmega(iterate,Omega), return; end
end

% Set up linear program and remove non-finite rows
c = [zeros(n,1); 1];
b = [iterate.x; -iterate.x; Omega.u; -Omega.l];
A = [ eye(n) ,  -ones(n,1); ...
        -eye(n) ,  -ones(n,1); ...
        Omega.A,  zeros(m,1); ...
        -Omega.A,  zeros(m,1)];
ind      = find(~isfinite(b));
b(ind)   = [];
A(ind,:) = [];

% Solve LP and assess if new iterate is feasible
[iterate.x,dist,err] = linprog(c,A,b,[],[],[],[],[],optimset('Display','off'));
iterate.x(end) = [];
if (err <= 0)
    error('mads:initialPoint',['MADS unable to find feasible initial point ', ...
            '(MakeFeasible).']);
end
return;

%*******************************************************************************
% FUNCTIONS CALLED BY MADS MAIN LOOP
%*******************************************************************************
%*******************************************************************************
% search:  Performs the SEARCH step of the MADS/GPS algorithm.
%          A finite number of points are evaluated according to the k-th Search
%          strategy set by the variable Search(k).Type. (Default: no search)
% ------------------------------------------------------------------------------
% Called by: mads
% Calls:     lhs, ga, recalSurrogate, optimizeSurrogate, evalPointSet, poll, 
%            standDirections, updateOmega, gridsamp (DACE)
% VARIABLES:
%  success         = flag indicating success or failure of Poll step
%  Filter          = structure of filter data
%    .F            =   indices of Cache points in the filter
%  Problem         = structure of optimization problem data
%    .nameCache    =   name of the base workspace Cache variable
%    .Omega        =   structure defining linear constraints
%      .plist      =     list of allowable discrete variable values
%    .isMVP        =   flag indicating an MVP problem
%    .adaptOmega   =   flag for updating an MVP's Omega parameters
%  Options         = structure of MADS parameters
%    .pollComplete =   turns on/off evaluation of ALL Poll points
%    .computeGrad  =   flag for computing any available gradients
%    .SurOptimizer =   string containing name of surrogate optimizer
%    .Search       =   vector of structures of Search data
%    .nSearches    =   number fo different Searches used during run
%    .pollStrategy =   string identifying selected Poll strategy
%    .pollComplete =   flag for performing a complete Poll step
%  RunData         = structure of run statistics and parameters
%    .scale        =   scale factors for scaling of mesh directions
%    .pcenter      =   current Poll center
%      .n          =     number of continuous variables
%    .delta        =   mesh size parameter
%    .nIter        =   iteration count
%  Cache           = structure containing of data on all past iterates
%    .iterate      =   vector of previously processed iterates
%    .size         =   number of iterates in Cache
%    .tol          =   tolerance for determining if a point is in Cache
%  full            = temporary storage of Options.pollComplete
%  grad            = temporary storage of Options.computeGrad
%  optimizer       = temporary storage of Options.SurOptimizer
%  Cache           = temporary storage of Cache.iterate
%  tol             = temporary storage of Cache.tol
%  Search          = structure of Search data 
%    .type         =   string indicating which search method to use
%    .nIter        =   number of iterations to perform Search
%    .nPoints      =   number of search points to be evaluated
%    .file         =   optional name of user-defined Search file
%    .dace         =   structure of DACE surrogate data
%      .reg        =     regression function handle
%      .corr       =     correlation function handle
%      .theta      =     correlation function parameter
%      .isotropic  =     flag indicating isotropic (all equal) theta values
%    .nw           =   structure of NW surrogate data
%      .kernel     =     name of NW kernel function
%      .local      =     flag indicating if surrogate is restricted to a region
%  surrogate       =   structure of surrogate data
%    .recalibrator =     surrogate recalibration function filename
%    .evaluatror   =     surrogate evaluation function filename
%  pcenter         = temporary poll center
%  D               = Mesh directions used in LHS
%  S               = Search points to be evaluated
%  N               = discrete neighbors of LHS points
%  q               = number of mesh points in each direction
%  LB,UB           = lower/upper bounds on continuous variables
%  meshPoints      = mesh points for coarse mesh point evaluations
%  param           = cell array of parameters output from user surrogate file.
%*******************************************************************************
function [success,Filter] = search(Problem,Options,RunData,Filter)

Cache = getappdata(0,Problem.nameCache);
success = 0;

% Shortcuts
full      = Options.pollComplete;
grad      = Options.computeGrad;
optimizer = Options.SurOptimizer;
tol       = Cache.tol;

% Include only iterates with the same categorical variable values (John Dunlap)
if (Problem.isMVP)
    pID = getpID(RunData.pcenter.p,Problem.Omega.plist);
    ind = find(Cache.pID == pID);
    Cache.iterate     = Cache.iterate(ind);
    Cache.isSurrPoint = Cache.isSurrPoint(ind);
end
ind   = find(Cache.isSurrPoint);  
Cache = Cache.iterate(ind);

% Identify which Search will be performed
if isempty(Options.Search), return, end
Search      = Options.Search(1);
Search.dace = Options.dace(1);
Search.nw   = Options.nw(1);
for k = 2:Options.nSearches
    if (RunData.nIter >= sum([Options.Search(1:k-1).nIter]))
        Search      = Options.Search(k);
        Search.dace = Options.dace(k);
        Search.nw   = Options.nw(k);
    end
end

% Change Surrogate to LHS if Cache is empty
if length(Cache) <= 1
    switch upper(Search.type)
        case {'DACE','NW','CUSTOMS'}
            Search.type = 'LHS';
            Search.nIter = 1;
            Search.nPoints = 2*RunData.pcenter.n;
            Search.file = '';
    end 
end

% Initialize variables for DACE, NW, or Custom surrogates
surrogate = [];
if isappdata(0,'SUR'), surrogate = getappdata(0,'SUR'); end

% Call appropriate Search functions
switch upper(Search.type)
    
    % Do nothing
    case {'NONE'}
        return
        
        % Poll around best n filter points
    case {'SPOLLI','CPOLLI','GPOLLI'}
        Options.pollStrategy = 'Standard_2n';
        if strcmp(Search.type,'CPollI')
            Options.pollComplete = 1;
        end
        if strcmp(Search.type,'GPollI') && Options.computeGrad
            Options.pollStrategy = 'Gradient_3n_L2';
        end
        for k = 1:min(length(Filter.F),Search.nPoints)
            pcenter = Cache(Filter.F(k));
            if (Problem.adaptOmega)
                [Problem.Omega,RunData.scale] = updateOmega(Problem,Options,pcenter);
            end
            [success,Filter] = poll(Problem,Options,RunData,pcenter,Filter);
            if (success), break, end
        end
        return
        
        % Latin Hypercube Search on the mesh
    case 'LHS'
        n = RunData.pcenter.n;
        Options.removeRedundancy = 0;
        [Problem.Omega,RunData.scale] = updateOmega(Problem,Options,RunData.pcenter);
        scale = min(RunData.scale,RunData.scale/norm(RunData.scale,inf));
        sdelta = RunData.delta * scale;
        LB = fix((Problem.Omega.l(1:n) - RunData.pcenter.x)./sdelta);
        UB = fix((Problem.Omega.u(1:n) - RunData.pcenter.x)./sdelta);
        maxrand = 2^5;
        LB(LB < -maxrand) = -maxrand/2;
        UB(UB >  maxrand) =  maxrand/2;
        z = lhs(UB-LB+1,Search.nPoints,2);
        S = [];
        if ~isempty(z)
            z = z + repmat(LB,1,size(z,2));
            for k = 1:size(z,2)
                S(k).x = RunData.pcenter.x + sdelta.*z(:,k);
                S(k).p = RunData.pcenter.p;
            end
        end
        full = (RunData.nIter == 0);
        
        % Deterministic Sampling on a Coarse Mesh (requires DACE package)
    case 'MESH'
        n  = RunData.pcenter.n;
        q  = ceil(Search.nPoints^(1/n));
        Options.removeRedundancy = 0;
        [Problem.Omega,RunData.scale] = updateOmega(Problem,Options,RunData.pcenter);
        LB = Problem.Omega.l(1:n)';
        UB = Problem.Omega.u(1:n)';
        LB(~isfinite(LB)) = -1e+16;
        UB(~isfinite(UB)) =  1e+16;
        meshPoints = gridsamp([LB;UB],q);
        for k = 1:size(meshPoints,1)
            S(k).x = meshPoints(k,:)';
            S(k).p = RunData.pcenter.p;
        end
        
        % Genetic Algorithm
    case 'GA'
        Options.removeRedundancy = 0;
        [Problem.Omega,RunData.scale] = updateOmega(Problem,Options,RunData.pcenter);
        S = ga(Problem,Options,Search.nPoints,RunData,Cache,Filter);
        
        % DACE Surrogate
    case 'DACE'
        surrogate.recalibrator = 'dacefit';
        surrogate.evaluator    = 'predictor';
        if Search.merit
            Search.merit = Search.merit/(2^RunData.nIter);
        end
        surrogate.merit = Search.merit;
        surrogate.local = Search.local;
        reg  = str2func(Search.dace.reg);
        corr = str2func(Search.dace.corr);
        if (Search.dace.isotropic)
            n = strcmp(Search.dace.corr,'correxpg') + 1;
        else
            n = strcmp(Search.dace.corr,'correxpg') + Cache(end).n;
        end
        theta = ones(n,1)*Search.dace.theta;
        surrogate = recalSurrogate(Cache,surrogate,optimizer,Filter,reg,corr,theta);
        S = optimizeSurrogate(Problem,RunData.pcenter,surrogate,optimizer, ...
            Search.nPoints,tol);
        
        % Nadaraya-Watson Surrogate
    case 'NW'
        surrogate.recalibrator = 'buildNW';
        surrogate.evaluator    = 'evalNW';
        if Search.merit
            Search.merit = Search.merit/(2^RunData.nIter);
        end
        surrogate.merit = Search.merit;
        surrogate.local = Search.local;
        surrogate = recalSurrogate(Cache,surrogate,optimizer,Filter, ...
            Search.nw.kernel,Search.nw.sigma, ...
            Search.nw.lower, Search.nw.upper);
        S = optimizeSurrogate(Problem,RunData.pcenter,surrogate,optimizer, ...
            Search.nPoints,tol);
        
        % Custom surrogate method
    case 'CUSTOMS'
        if (~exist(Search.file,'file'))
            error('mads:search:file',[Search.file '.m not found (search).']);
        end
        if upper(strcmp(optimizer,'CUSTOM'))
            S = feval(Search.file,Problem,Cache,Search,surrogate,tol);
        else
            if Search.merit
                Search.merit = Search.merit/(2^RunData.nIter);
            end
            surrogate.merit   = Search.merit;
            surrogate.local   = Search.local;
            [surrogate,param] = feval(Search.file);
            surrogate = recalSurrogate(Cache,surrogate,optimizer,Filter,param{:});
            S = optimizeSurrogate(Problem,RunData.pcenter,surrogate,optimizer, ...
                Search.nPoints,tol);
        end
        
        % Custom Search method
    case 'CUSTOM'
        if (~exist(Search.file,'file'))
            error('mads:search:file',[Search.file '.m not found (search).']);
        end
        S = feval(Search.file,Problem,Options,RunData.delta,RunData.pcenter);
        
        % Invalid Search option
    otherwise
        S = [];
        error('mads:search:choice','Invalid choice of Search Method (search).');
end

[S,success,Filter,RunData.RS] = evalPointSet(Problem,S,full,Options,Filter, ...
    RunData.RS,2);
return;

%*******************************************************************************
% poll: Performs the POLL step of the MADS algorithm.  Given a Poll center, a
%       set of Poll direction vectors, and a mesh size parameter, the
%       neighboring mesh points are evaluated, seeking an improved mesh point.
% ------------------------------------------------------------------------------
% Called by: mads, Search, mvpPoll
% Calls:     standDirections, gradDirections, getPollOrder, evalPointSet
% VARIABLES:
%  success         = flag indicating success or failure of Poll step
%  Filter          = structure describing the filter
%  Problem         = structure describing the optimization problem
%    .Omega        =   structure defining Omega = {x: l <= A*x <= u}
%    .adaptOmega   =   flag for updating an MVP's Omega parameters
%  Options         = structure containing various MADS parameters
%    .computeGrad  =   flag for computing any available gradients
%    .pollStrategy =   string identifying selected Poll strategy
%    .pollOrder    =   string identifying selected Poll order strategy
%    .pollComplete =   flag for performing a complete Poll step
%    .tolBind      =   active constraint tolerance
%  RunData         = structure of run statistics and parameters
%    .scale        =   scale factors for scaling mesh directions
%    .delta        =   mesh size parameter
%    .goodD        =   index of successful Poll direction
%    .porder       =   current poll order
%  pcenter         = center of Poll set, around which Poll is performed
%    .x            =   continuous variables
%    .p            =   categorical variables
%    .h            =   constraint violation function value
%    .gradf        =   gradient of f
%    .gradh        =   gradient of h
%  D               = matrix whose columns form the direction vectors
%  g               = gradient direction with which to prune
%  infGrad         = indices of non-finite partial derivatives
%  donef           = indices of zero partial derivatives of f
%  doneh           = indices of zero partial derivatives of h
%  surrogate       = structure of surrogate model data
%  order           = order in which directions will be polled
%  P               = The Poll set
%*******************************************************************************
function [success,Filter,varargout] = poll(Problem,Options,RunData,pcenter, ...
    Filter)

% A gradient poll without gradients available becomes a standard poll
if (Options.computeGrad && ~isfield(pcenter, 'gradf'))
    Options.pollStrategy = 'Standard_2n';
end

% Retrieve scaled Poll directions
switch Options.pollStrategy
    case {'Standard_2n','Standard_n+1','Custom_2n','Custom_n+1'}
        D = standDirections(pcenter.x,Problem.Omega,Options.pollStrategy, ...
            Options.tolBind,RunData.scale);
    case {'MADS_2n','MADS_n+1'}
        D = madsDirections(pcenter.x,Options.pollStrategy,RunData.delta, ...
            RunData.scale);
    case {'Gradient_2n','Gradient_n+1'}
        D = standDirections(pcenter.x,Problem.Omega,Options.pollStrategy, ...
            Options.tolBind,RunData.scale);
        if (pcenter.h == 0)
            g = pcenter.gradf;
        else
            g = pcenter.gradh;
        end
        infGrad = ~isfinite(g);
        g(infGrad) = 0;
        D(:,(infGrad'*D == 0 & g'*D > 0.0)) = [];
    case {'Gradient_3n_L1','Gradient_3n_L2','Gradient_3n_LInf','Gradient_3n2n'}
        D = [];
        donef = all(pcenter.gradf == 0.0);
        doneh = (isempty(pcenter.gradh) || all(pcenter.gradh) == 0.0);
        if ~(donef && doneh)
            D = gradDirections(pcenter,Problem.Omega,Options.pollStrategy,Filter, ...
                RunData.delta,Options.tolBind,RunData.scale);
        end
    otherwise
        error('mads:poll:choice','Invalid Poll strategy type (poll).');
end

% Construct the POLL set
P = [];
for k = 1:size(D,2)
    P(k).x = pcenter.x + RunData.delta*D(:,k);
    P(k).p = pcenter.p;
end

% Get surrogate information for setting surrogate Poll order
surrogate = [];
if isappdata(0,'SUR')
    surrogate = getappdata(0,'SUR');
end
if strcmp(Options.pollOrder,'Surrogate') && isempty(surrogate)
    Options.pollOrder = 'DynamicRanked';
end

% Set the polling order, and evaluate the POLL points
order = getPollOrder(Options.pollOrder,D,P,RunData,surrogate);
P = P(order);

Problem.adaptOmega = 0;
[P,success,Filter,RunData.RS] = evalPointSet(Problem,P,Options.pollComplete, ...
    Options,Filter,RunData.RS,1);

% Update parameters need for dynamic poll ordering
if (success)
    RunData.goodD  = D(:,success);
    RunData.porder = order;
    success = 1;
end
varargout = {RunData};
return;

%*******************************************************************************
% mvpPoll:  Performs discrete neighbor and extended polls for MVP problems.
% ------------------------------------------------------------------------------
% Called by: mads
% Calls:     evalPointSet, poll, getPollCenter, createFilter, updateOmega
% VARIABLES:
%  success         = flag indicating success or failure of MVP Poll step
%  Problem         = structure describing the optimization problem
%    .nameCache    =   name of the base workspace Cache variable
%    .File.N       =   name of discrete neighbor file
%    .Omega        =   structure defining linear constraints
%      .plist      =     list of allowable discrete variable values
%    .adaptOmega   =   flag for updating an MVP's Omega parameters
%  Options         = structure containing various MADS parameters
%    .pollComplete =   flag for performing a complete Poll step
%    .pollCenter   =   code for identifying the Poll center
%  RunData         = structure containing MADS run statistics
%    .delta        =   current mesh size
%    .pcenter      =   current Poll center
%    .fxi          =   f-value threshold for triggering Extended Poll
%    .scale        =   scale factors for scaling of mesh directions
%  Cache           = database of all previously computed iterates
%    .Filter       =   structure of parameters describing the main filter
%      .hmin       =     minimum h-value of an infeasible point
%      .hmax       =     maximum allowable h-value
%  N               = discrete neighbors of the current poll center
%    .f            =   objective function value
%    .h            =   constraint violation function value
%  BestF           = best feasible iterate
%    .f            =   objective function value
%  BestI           = least infeasible iterate
%    .h            =   constraint violation function value
%  fplusxi         = BestF.f + RunData.fxi
%  hplusxi         = BestF.h + RunData.hxi (or hmax if this is too large)
%  ePollF          = indices of N triggering Extended Poll due to f
%  ePollH          = indices of N triggering Extended Poll due to h
%  ePoll           = indices of N triggering Extended Poll due to f or h
%  new             = temporary storage
%  order           = sorted index for N by f-value or h-value
%  pcenter         = current Extended Poll center
%  nx              = number of continuous variables in pcenter
%  Filter          = structure describing the temporary local filter
%  ePsuccess       = flag indicating a successful Extended Poll step
%  unfiltered      = flag indicating an unfiltered point has been found
%*******************************************************************************
function success = mvpPoll(Problem,Options,RunData)

% Retrieve set of discrete neighbors
N = feval(Problem.File.N,Problem,RunData.pcenter,Problem.Omega.plist, ...
    RunData.delta);
if isempty(N)
    success = 0;
    return;
end

% Evaluate set of discrete neighbors
Cache = getappdata(0,Problem.nameCache);
[N,success,TempFilter,RunData.RS] = evalPointSet(Problem,N,1,Options, ...
    Cache.Filter,RunData.RS,1);
Cache = getappdata(0,Problem.nameCache);
Cache.Filter = TempFilter;
setappdata(0,Problem.nameCache,Cache);
N = N(find(~isempty([N.f])));
success = ~~success;

% If unsuccessful, begin extended poll around "good" discrete neighbors
if (~success)
    
    % Determine candidate poll centers for extended polling
    BestF = getPollCenter(0,Cache.Filter,Cache);
    BestI = getPollCenter(1,Cache.Filter,Cache);
    fplusxi = BestF.f + RunData.fxi;
    hplusxi = min(BestI.h + RunData.hxi,Cache.Filter.hmax);
    ePollF = find([N.h] <  Options.hmin & [N.f] < fplusxi);
    ePollH = find([N.h] >= Options.hmin & [N.h] < hplusxi);
    
    % Determine the order in which extended poll centers get polled around
    if (~isempty(ePollF))
        [new,order] = sort(N(ePollF).f);
        ePollF = ePollF(order);
    end
    if (~isempty(ePollH))
        [new,order] = sort(N(ePollH).h);
        ePollH = ePollH(order);
    end
    if (RunData.pcenter.h < Cache.Filter.hmin)
        ePoll = [ePollF, ePollH];
    else
        ePoll = [ePollH, ePollF];
    end
    
    % Search/Poll around each selected discrete neighbor using local filter
    for k = 1:length(ePoll)
        
        % Create MVP filter and populate it with the discrete neighbor
        pcenter = N(ePoll(k));
        Filter = createFilter(Options.hmin,hplusxi,1,0,[]);
        [unfiltered,Filter] = updateFilter(pcenter,Filter,Cache);      
        if (Problem.adaptOmega)
            [Problem.Omega,RunData.scale] = updateOmega(Problem,Options,pcenter);
        end
        
        % Begin EXTENDED SEARCH step using same Search type is in the SEARCH step
        TempRunData = RunData;
        TempRunData.pcenter = pcenter;
        TempOptions = Options;
        for j = 1:Options.nSearches
            TempOptions.Search(j) = Options.Search(j);
        end
        Cache        = getappdata(0,Problem.nameCache);
        startCache   = Cache.size + 1;
        [success,TempFilter] = search(Problem,TempOptions,TempRunData, ...
            Cache.Filter);
        Cache        = getappdata(0,Problem.nameCache);
        endCache     = Cache.size;
        Cache.Filter = TempFilter;
        setappdata(0,Problem.nameCache,Cache);
        
        % If EXTENDED SEARCH was unsuccessful, test versus local filter      
        if (~success)
            for j = startCache:endCache
                [unfiltered,Filter] = updateFilter(Cache.iterate(j),Filter,Cache);
            end
            pcenter = getPollCenter(Options.pollCenter,Filter,Cache);
        end
        
        % Begin EXTENDED POLL
        while (~success)
            startCache = Cache.size + 1;
            [ePsuccess,Filter] = poll(Problem,Options,RunData,pcenter,Filter);
            Cache = getappdata(0,Problem.nameCache);
            endCache = Cache.size;
            
            % If Extended Poll was unsuccessful, try next discrete neighbor
            if (~ePsuccess), break; end
            
            % If Extended Poll step was successful, determine new incumbent
            if (Options.pollComplete)
                for j = startCache:endCache
                    iterate = Cache.iterate(j);
                    [unfiltered,TempFilter] = updateFilter(iterate,Cache.Filter,Cache);
                    if (unfiltered)
                        success = j;
                        Cache.isSurrPoint(j) = 1;
                    end
                    Cache.Filter = TempFilter;
                    setappdata(0,Problem.nameCache,Cache);
                end
            else
                [iterate,iCache] = getPollCenter(Options.pollCenter,Filter,Cache);
                if Options.runStochastic, Cache.size = iCache; end
                [unfiltered,TempFilter] = updateFilter(iterate,Cache.Filter,Cache);
                if (unfiltered)
                    success = k;
                    Cache.isSurrPoint(iCache) = 1;
                end
                Cache.Filter = TempFilter;
                setappdata(0,Problem.nameCache,Cache);
            end
        end
        if (success && ~Options.pollComplete), break; end
    end
end
return

%*******************************************************************************
% update:  Updates MADS run parameters.
% ------------------------------------------------------------------------------
% Called by: mads, processInput
% Calls:     getPollCenter, updateOmega, plotHistory
% VARIABLES:
%  Problem             = structure defining optimization problem
%    .nameCache        =   name of the base workspace Cache variable
%    .Omega            =   structure defining linear constraints
%    .isMVP            =   flag indicating an MVP problem
%    .adaptOmega       =   flag for updating an MVP's Omega parameters
%  RunData             = structure holding MADS run statistics
%    .delta            =   current mesh size
%    .nIter            =   iteration count
%    .nFunc            =   cumulative number of function evaluations
%    .grad             =   cumulative number of gradient evaluations
%    .time             =   cumulative CPU time used
%    .nFails           =   current count of consecutive Poll failures
%    .nCacheHits       =   cumulative number of Cache hits
%    .nFunc0           =   function evaluations during initialization
%    .pcenter          =   current poll center
%    .fxi              =   f-value threshold for triggering Extended Poll
%    .hxi              =   h-value threshold for triggering Extended Poll
%    .goodD            =   index of the most recent successful direction
%    .stopRun          =   flag for stopping run immediately
%    .scale            =   scale factors for each continuous variable
%  Options             = structure holding MADS parameters
%    .computeGrad      =   flag for computing any available gradients
%    .countCache    =   flag for counting Cache points as function calls
%    .delta0           =   initial mesh size
%    .deltaMax         =   maximum mesh size 
%    .meshRefine       =   mesh size refinement factor
%    .meshCoarsen      =   mesh size coarsening factor
%    .accelerate       =   flag for decreasing the mesh refinement factor
%    .ePollTriggerF    =   f-value trigger for executing extended POLL
%    .ePollTriggerH    =   h-value trigger for executing extended POLL
%    .pollCenter       =   code indicating which poll center to select
%    .plotHistory2     =   turns on/off plot of f-value vs. #f-evals
%    .plotColor        =   string indicating color of history plot line
%    .stophandle       =   handle for external object that terminates run
%    .runUntilFeasible =   flag for running MADS only until feasible
%    .runOneIteration  =   flag for running MADS one iteration at a time
%  success             = flag indicating iteration successful/failure
%  Cache               = database of all previously computed iterates
%    .size             =   number of iterates currently in the Cache
%    .nHits            =   number of Cache hits
%  Filter              = structure containing filter data
%    .hmax             =   maximum allowable constraint violation h-value
%    .feasible         =   indices of best feasible points in the Cache
%  meshScaleFactor     = factor for increasing or decreasing mesh size
%  BestF               = best feasible iterate found thus far
%*******************************************************************************
function [Problem,RunData] = update(Problem,Options,RunData,success)

Cache = getappdata(0,Problem.nameCache);

% Update run statistics for either initial or non-initial iterates
if (success < 0)
    RunData.delta       = Options.delta0;
    RunData.nIter       = 0;
    RunData.nFunc0      = (~Options.countCache)*Cache.size;
    RunData.nFunc       = 0;
    RunData.grad        = Options.computeGrad;
    RunData.time        = 0;
    RunData.nFails      = 0;
    RunData.meshRefine  = Options.meshRefine;
    RunData.meshCoarsen = Options.meshCoarsen;
    RunData.hxi         = Options.ePollTriggerH*Cache.Filter.hmax;
else
    meshScaleFactor    = ~~success*RunData.meshCoarsen + ...
        ~success*RunData.meshRefine;
    RunData.delta      = min(meshScaleFactor*RunData.delta, Options.deltaMax);
    RunData.nIter      = RunData.nIter + 1;
    RunData.nFunc      = sum(Cache.nFunc(1:Cache.size)) - RunData.nFunc0;
    RunData.grad       = RunData.grad + (Options.computeGrad && success);
    RunData.time       = toc;
    RunData.nFails     = ~success*(RunData.nFails + 1);
    RunData.meshRefine = RunData.meshRefine/(1+(~success && Options.accelerate));
end
if Options.runStochastic
    RunData.RS.nFuncLeft = Options.Term.nFunc - RunData.nFunc;
end
RunData.pcenter      = getPollCenter(Options.pollCenter,Cache.Filter,Cache);
RunData.nCacheHits   = Cache.nHits;

% Update StopRun flag, as appropriate
RunData.stopRun = 0;
if (~isempty(Options.stophandle))
    RunData.stopRun = get(Options.stophandle,'UserData');
end
if (Options.runUntilFeasible && ~isempty(Cache.Filter.feasible))
    RunData.stopRun = 1;
end
if (Options.runOneIteration)
    RunData.stopRun = 1;
end
drawnow

% For MVP, update Omega, scale factors, and extended poll thresholds
if (Problem.isMVP && success)
    if (Problem.adaptOmega)
        [Problem.Omega,RunData.scale] = updateOmega(Problem,Options, ...
            RunData.pcenter);
    end
    BestF = getPollCenter(0,Cache.Filter,Cache);
    RunData.fxi = Options.ePollTriggerF*max(Options.Term.delta,abs(BestF.f));
end
if (Options.plotHistory2)
    plotHistory(Options.hplothandle,Options.plotColor,Cache);
end
return

%*******************************************************************************
% FUNCTIONS CALLED BY MADS SEARCH ROUTINE
%*******************************************************************************
%*******************************************************************************
% lhs:  Perform a Latin Hypercube Search for integer vectors
% ------------------------------------------------------------------------------
% Called by: search
% VARIABLES:
%  z        = matrix whose columns are LH vectors
%  ngrid    = vector whose elements are the number of points in each dimension
%  budget   = the maximum total number of points to be sampled
%  strength = stength (desired number of points per dimension) of the LHS
%  nPoints  = actual number of points generated
%  nleft    = number of indices in each dimension not already chosen
%  n        = dimension of sample space
%  ilist    = indices of already chosen numbers in each direction
%  endflag  = flag indicating usage of all points in a direction
%*******************************************************************************
function [z] = lhs(ngrid,budget,strength)

% Initialize variables and determine the number of points to be evaluated
n       = length(ngrid);
[z,ilist{1:n}] = deal([]);
nPoints = min(strength*max(ngrid), prod(ngrid));
nPoints = min(nPoints,budget);
nleft   = ngrid;

% Begin main loop to generate LHS points
for i = 1:nPoints
    
    % Generate random vector of integers between 0 and nleft
    z(:,i) = fix(nleft.*rand(n,1) - eps);
    
    % Adjust and update indices, referencing only unused ones
    for k = 1:n
        endflag = 1;
        for j = 1:length(ilist{k})
            if (z(k,i) >= ilist{k}(j))
                z(k,i) = z(k,i) + 1;
            else
                endflag = 0;            
                break;
            end
        end
        if isempty(j), j = 0; end
        ilist{k} = [ilist{k}(1:j-1+endflag), z(k,i), ilist{k}(j+endflag:end)];
        if (length(ilist{k}) >= ngrid(k)), ilist{k} = []; end
        nleft(k) = ngrid(k) - length(ilist{k});
    end
end
return;

%*******************************************************************************
% ga:  Perform Genetic Algorithm on a mesh, as part of Search.
% ------------------------------------------------------------------------------
% Called by: search
% Calls:     setPenaltyTerm
% VARIABLES:
%  S         = Set of "good" iterates to be evaluated
%  Problem   = structure of optimization problem data
%  Options   = structure of user options
%  nPoints   = number of points to return
%  pcenter   = current Poll center
%  iterate   = vector of iterates
%  Filter    = index into iterate of filter and best feasible points
%  arg       = temporary string used in writing penalty function file
%  n         = number of continuous variables
%  a         = parameter used in penalty function used in surrogate
%  f         = name of penalty function
%  fid       = file ID (handle) of penalty function
%  LB,UB     = storage of upper and lower variable bounds for GA algorithm use
%  options   = options for GA optimizer
%  es,result = GA optimizer output variables
%*******************************************************************************
function S = ga(Problem,Options,nPoints,RunData,iterate,Filter)

% Write penalty function file
if Problem.isMVP
    setappdata(0,'P',RunData.pcenter.p);
    arg = '(x,p);';
else
    arg = '(x);';
end
a   = setPenaltyTerm(iterate,Filter);
n   = RunData.pcenter.n;
f   = 'gaPenalty';
fid = fopen([f,'.m'], 'w');
fprintf(fid,'%s\n',    ['function [z] = ', f, '(x);'           ], ...
    ['n = ', int2str(n),';'                 ], ...
    'if isappdata(0,''P'')'                 , ...
    '   p = getappdata(0,''P'');'           , ...
    'end'                                   , ...
    ['[fx,cx] = ', Problem.File.F, arg      ]);
if exist(Problem.File.O,'file')
    fprintf(fid,'%s\n', ['[A,l,u] = ', Problem.File.O, '(n);'   ], ...
        'A = [A(n+1:end,:); -A(n+1:end,:)];'    , ...
        'b = [u(n+1:end,:); -l(n+1:end,:)];'    , ...
        'ind = ~isfinite(b);'                   , ...
        'b(ind)   = [];'                        , ...
        'A(ind,:) = [];'                        , ...
        'cx = [cx ; A*x - b];'                  );
end
fprintf(fid,'%s\n',     'cplus   = (cx > 0).*cx;'               , ...
    'hx      = norm(cplus)^2;'              , ...
    ['z       = fx + ', num2str(a), '*hx;'  ], ...
    'return'                               );
fclose(fid);
rehash;

% Set up bounds for the GA code (requires full set of finite variable bounds)
LB = Problem.Omega.l(1:n);
UB = Problem.Omega.u(1:n);
LB(~isfinite(LB)) = -1/eps;
UB(~isfinite(UB)) =  1/eps;

% Call GA code
options = es_options_new('TolX',Options.Term.delta,'SigmaFacStart', 1, ...
    'LBound',LB,'UBound',UB,'scaling',RunData.scale);
es      = es_new(n, LB, UB, options);
es      = es_run(es, f, nPoints);               % minimize f with an ES struct
result  = es_get(es, 'result');                 % get result
S.x = result.x;                                 % return current best point
S.p = RunData.pcenter.p;
delete([f,'.m']);

return

%*******************************************************************************
% setPenaltyTerm: Get penalty parameter for applying GA to constrained problem.
% ------------------------------------------------------------------------------
% Called by: ga, recalSurrogate
% VARIABLES:
%  a           = parameter used in penalty function used in surrogate or GA
%  iterate     = vector of iterates
%    .f        =   objective function value
%    .h        =   constraint violation function value
%  Filter      = index into iterate of filter and best feasible points
%    .F        =   indices of filter points
%    .feasible =   indices of best feasible solutions
%  FF          = objective function values of filter points
%  HH          = constraint violation function values of filter points
%*******************************************************************************
function a = setPenaltyTerm(iterate,Filter)

% Compute penalty parameter
FF = [iterate(Filter.F).f]';
HH = [iterate(Filter.F).h]';
if isempty(Filter.feasible)
    a = 1000;
else
    FF = [iterate(Filter.feasible(1)).f; FF]';
    HH = [iterate(Filter.feasible(1)).h; HH]';
    a  = max((FF(1:end-1)-FF(2:end))./(HH(2:end)-HH(1:end-1)));
    if isempty(a), a = 0.5; end
end
return

%*******************************************************************************
% recalSurrogate:  Construct or recalibrate a surrogate.
% ------------------------------------------------------------------------------
% Called by: search
% Calls:     setPenaltyTerm
% Calls: dacefit (DACE), buildNW (NW), or other surrogate-specific recalibrator
% VARIABLES:
%  iterate     = vector of iterates (all in the Cache)
%  surrogate   = structure containing surrogate information
%  optimizer   = string identifying the optimizer used to find S
%  Filter      = indices into iterate that contain the filter
%  X           = continuous variables values of iterates used for surrogate
%  F,C         = objective and constraint values of iterates used for surrogate
%  H           = constraint violation function values of surrogate iterates
%  n           = number of points used to build surrogate
%  nc          = number of nonlinear constraints
%  FF,HH       = objective and constraint values of filter points
%  a           = parameter used in penalty function used in surrogate
%  distfactor  = factor used in defining a distance parameter
%  maxDistance = largest distance between any two data sites
%*******************************************************************************
function surrogate = recalSurrogate(iterate,surrogate,optimizer,Filter,varargin)

% Set up data sites and recalibrate surrogate
X  = [iterate.x]';
F  = [iterate.f]';
H  = [iterate.h]';
n  = length(iterate);
nc = length(iterate(1).c);
C  = zeros(n,nc);
for k = 1:n
    if ~isempty(iterate(k).c)
        C(k,:) = iterate(k).c;
    end
end
if strcmp(optimizer,'GA')
    a = setPenaltyTerm(iterate,Filter);
    [surrogate.f] = feval(surrogate.recalibrator,X,F+a*H,varargin{:});
else
    [surrogate.f] = feval(surrogate.recalibrator,X,F,varargin{:});
    if (length(iterate(end).c) > 0)
        surrogate.c = feval(surrogate.recalibrator,X,C,varargin{:});
    end
end

% Compute trust region radius and distance and merit function penalty term
if surrogate.local
    maxDistance = 0;
    for i = 1:length(F);
        for j = 1:i-1
            distance = norm(X(i,:) - X(j,:));
            if (distance > maxDistance)
                maxDistance = distance;
            end
        end
    end
else
    maxDistance = inf;
end
surrogate.trust = maxDistance/2;
surrogate.dist  = 10*surrogate.merit*(max(F) - min(F));
surrogate.X = X;
surrogate.F = F;
return

%*******************************************************************************
% optimizeSurrogate:  Optimize a surrogate problem.
% ------------------------------------------------------------------------------
% Called by: search
% Calls: dacefit (DACE), buildNW (NW), fmincon (Optimization Toolbox)
% VARIABLES:
%  S          = set of points generated to be evaluated
%  Problem    = structure containing optimation problem data
%    .Omega   =   structure defining feasible region
%      .A     =     matrix of linear coefficients
%      .l     =     vector of lower bounds
%      .u     =     vector of upper bounds
%  pcenter    = current Poll center
%  surrogate  = structure containing surrogate information
%  optimizer  = string identifying the optimizer used to find S
%  nPoints    = number of points to be returned by the optimizer
%  tol        = desired accuracy of optimizer
%  f          = name of the surrogate objective function file
%  C          = name of the surrogate constraints function file
%  fid        = file ID (handle) for a surrogate problem file
%  Sur        = structure of input data for running MADS on the surrogate
%  SurRunData = RunData for MADS run on the surrogate problem
%  SurCache   = Cache for the MADS run on the surrogate problem
%  temp1,2    = temporary storage
%  n          = number of continuous variables
%  options    = options for GA optimizer
%  es,result  = GA optimizer output variables
%  sf         = surrogate function values of points from GA optimizer
%  order      = indices of GA-produced iterates, sorted by sf
%  y          = temporary storage
%*******************************************************************************
function S = optimizeSurrogate(Problem,pcenter,surrogate,optimizer,nPoints,tol)

setappdata(0,'SUR',surrogate);

% Construct and store MATLAB commands for computing a surrogate function value
p{1}     = 'surrogate = getappdata(0,''SUR'');';
p{end+1} = 'isLocal = 0;';
p{end+1} = 'minDist = Inf;';
p{end+1} = 'for k = 1:size(surrogate.X,1)';
p{end+1} = '   normterm = norm(transpose(x) - surrogate.X(k,:));';
p{end+1} = '   minDist  = min(minDist, normterm);';
p{end+1} = '   if (normterm <= surrogate.trust)';
p{end+1} = '      isLocal = 1;';
p{end+1} = '      break';
p{end+1} = '   end';
p{end+1} = 'end';
p{end+1} = 'fx = 1/eps;';
p{end+1} = 'if isLocal';
p{end+1} = '   fx = feval(surrogate.evaluator,x,surrogate.f);';
p{end+1} = '   fx = fx - minDist*surrogate.dist;';
p{end+1} = 'end';
p{end+1} = 'if ~isfinite(fx), fx = 1/eps; end';

switch optimizer
    
    % Call optimizer FMINCON    
    case 'fmincon'
        
        % Write surrogate function files if they do not already exist
        f   = 'SurObj';
        if ~exist([f,'.m'],'file') 
            fid = fopen([f,'.m'], 'w');
            fprintf(fid,'%s\n',['function [fx] = ', f, '(x);'                  ]);
            fprintf(fid,'%s\n',p{:});
            fprintf(fid,'%s\n','return');
            fclose(fid);
            rehash
        end
        C = [];
        if isfield(surrogate,'c')
            C   = 'SurNLConst';
            if ~exist([C,'.m'],'file')
                fid = fopen([C,'.m'],'w');
                fprintf(fid,'%s\n', ...
                    ['function [C,Ceq] = ', C, '(x);'                 ], ...
                    'surrogate = getappdata(0,''SUR'');'              , ...
                    'if isfield(surrogate,''c'')'                     , ...
                    '   C = feval(surrogate.evaluator,x,surrogate.c);', ...
                    'end'                                             , ...
                    'Ceq = [];'                                       , ...
                    'return'                                         );
                fclose(fid);
                rehash
            end
        end
        
        % Call fmincon optimizer
        warning('off','all');
        if isempty(Problem.Omega.A) && isempty(C)
            S.x = fminunc(f,pcenter.x,optimset('Display','off','TolX',tol));
        else
            S.x = fmincon(f,pcenter.x,[Problem.Omega.A; -Problem.Omega.A], ...
                [Problem.Omega.u; -Problem.Omega.l],[],[],[],[], ...
                C,optimset('Display','off','TolX',tol));
        end
        warning('on','all');
        S.p = pcenter.p;
        
        % Call optimizer FMINCON
    case 'mads'
        
        % Setup options for Surrogate optimization by MADS
        Sur.Defaults          = mads_defaults('Surrogate');
        Sur.Problem           = Problem;
        Sur.Problem.nameCache = Sur.Defaults.nameCache;
        Sur.Problem.File.F    = [Problem.File.F, '_Sur'];
        Sur.Problem.File.C    = [Problem.File.C, '_Sur_Cache'];
        
        % Write surrogate function file for MADS if it does not already exist
        f   = Sur.Problem.File.F;
        if ~exist([f,'.m'],'file')
            fid = fopen([f,'.m'], 'w');
            fprintf(fid,'%s\n', ['function [fx,cx] = ', f, '(x);']);
            fprintf(fid,'%s\n',p{:});
            if isfield(surrogate,'c')
                fprintf(fid,'%s\n','cx = feval(surrogate.evaluator,x,surrogate.c);');
            else
                fprintf(fid,'%s\n','cx = [];');
            end
            fprintf(fid,'%s\n',   'return');
            fclose(fid);
            rehash
        end
        
        % Pass function to MADS optimizer
        [temp1,temp2,SurRunData,SurCache] = mads(Sur.Problem,pcenter, ...
            Sur.Defaults.Options);
        
        % Retrieve and evaluate a number of points
        S = [SurCache.iterate(SurCache.Filter.feasible(1:nPoints)), ...
                SurCache.iterate(SurCache.Filter.F(1:nPoints))];
        
        % Call Genetic Algorithm package as the surrogate optimizer
    case 'cmaes'
        LB = Problem.Omega.l(1:pcenter.n);
        UB = Problem.Omega.u(1:pcenter.n);
        options = es_options_new('TolX',tol,'SigmaFacStart',1,...
            'LBound',LB,'UBound',UB);
        sf = zeros(nPoints,1);
        for k = 1:nPoints
            es = es_new(pcenter.n, LB, UB, options);
            es = es_run_mod(es, surrogate.evaluator, 100, {surrogate.f}); %minimize
            result = es_get(es, 'result');
            S(k).x = result.x;
            S(k).p = pcenter.p;
            sf(k)  = feval(surrogate.evaluator,S(k).x,surrogate.f);
        end
        [y,order] = sort(sf);              % Sort points by surrogate f-value
        S = S(order);
end
return

%*******************************************************************************
% FUNCTIONS CALLED BY MADS POLL ROUTINE
%*******************************************************************************
%*******************************************************************************
% standDirections: Returns directions that positively span the tangent cone
%    at the current iterate, with respect to bound and linear constraints.
% ------------------------------------------------------------------------------
% Called by: search, poll
% Calls:     activeConstraints, getTangentCone, scaleDirections
% VARIABLES:
%  D        = matrix whose columns form the positive spanning set
%  x        = continuous variable values of the current iterate
%  Omega    = structure whose fields define Omega = {x: l <= Ax <= u}
%    .A     =   matrix of coefficients for the linear constraints
%    .l     =   vector of lower bounds
%    .u     =   vector of upper bounds
%  strategy = code indicating chosen poll strategy
%  tol      = tolerance for considering constraints active
%  scale    = scale factors for scaling the appropriate Poll directions
%  Ax       = A*x
%  n        = number of columns of A
%  I        = identity matrix
%  W        = underlying basis of Poll directions
%  B,N      = matrices whose columns together span the tangent cone
%*******************************************************************************
function D = standDirections(x,Omega,strategy,tol,scale)

% Initialize variables
Ax = Omega.A*x;
n  = length(x);

% Get user-provided Poll basis, if provided
I = eye(n);
W = I;
if isappdata(0,'PARAM') && strncmp(strategy,'Custom',6)
    if isfield(Param,'PollBasis')
        Param = getappdata(0,'PARAM');
        W = Param.pollBasis;
    end
end
if ~((rank(W) == n) && all(size(W) == [n,n]))
    error('mads:user:dim', ...
        ['User Poll basis has incompatible size or is not a basis ',...
            '(standDirections).']);
end

% Compute B and N, and scale each column in N
if (isequal(Omega.A,I))   % Bound constraints
    lactive = activeConstraints( Ax, Omega.l, Omega.u,tol);
    uactive = activeConstraints(-Ax,-Omega.u,-Omega.l,tol);
    active  = find(lactive | uactive);
    B = I(:,active);
    N = I(:,setdiff(1:n,active));
    if isempty(active), N = W; end
    B = scaleDirections(B,scale);
    N = scaleDirections(N,scale);
else                      % General linear constraints
    [B,N] = getTangentCone(x,Omega,tol,scale);
end

% Form directions that positively span the tangent cone at x
switch strategy
    case {'Standard_n+1','Gradient_n+1','Custom_n+1'}
        D = [-sum(N,2) N  B -B];
    case {'Standard_2n', 'Gradient_2n', 'Custom_2n'}
        D = [N -N B -B];
    otherwise
        error('mads:poll:choice','Invalid Poll strategy type (standDirections).');
end
D(:,~any(D)) = [];
return;

%*******************************************************************************
% madsDirections:  Returns a random set of directions that positively span
%    the tangent cone, with respect to bound and linear constraints.
% ------------------------------------------------------------------------------
% Called by: poll
% Calls:     scaleDirections
% VARIABLES:
%  D        = matrix whose columns form the positive spanning set
%  x        = continuous variable values of the current iterate
%  strategy = code indicating chosen poll strategy
%  delta    = mesh size parameter
%  scale    = scale factors for scaling the appropriate Poll directions
%  n        = number of continuous variables
%  z        = reciprocal of the square root of the mesh size
%  Z        = a random permutation of the integers 1-n
%  L,M      = matrices used in forming the MADS directions
%  B        = basis of directions
%  V        = random permutation of the rows and columns of B
%*******************************************************************************
function D = madsDirections(x,strategy,delta,scale)

n = length(x);
z = 1/sqrt(delta);

% Construct lower triangular nonsingular basis matrix
Z = rand(n);
L = tril(fix((2*z+1)*Z - eps), -1);
M = diag(z*sign(diag(Z-0.5)));
M(~M) = 1;
B = L + M;

% Randomly permute rows and columns of B
V = B(randperm(n),randperm(n));

% Construct positive basis
switch strategy
    case 'MADS_2n'
        D = [V -V];
    case 'MADS_n+1'
        D = [V -sum(V,2)];
    otherwise
        error('mads:poll:choice','Invalid Poll strategy type (madsDirections).');
end
D = scaleDirections(D,scale);
return

%*******************************************************************************
% gradDirections:  Compute Poll directions for a gradient-pruned Poll.
% ------------------------------------------------------------------------------
% Called by: poll
% Calls:     activeConstraints, getTangentCone, dVector, scaleDirections
% VARIABLES:
%  D                 = set of directions generated by gradient-pruned polling
%  iterate           = current iterate
%  Omega             = structure describing feasible region
%  strategy          = gradient poll type
%  Filter            = structure of Cache indices of iterates in the filter
%    .hmax           = max allowable constraint violation of any filter point
%  delta             = mesh size parameter
%  tol               = tolerence for determining if a constraint is active
%  scale             = direction scale factors
%  p                 = identifier for type of gradient-based direction
%  infGrad           = flags for identifying unavailable partial derivatives
%  g                 = gradient vector
%  d                 = gradient-based descent vector
%  B                 = indices of numerically binding nonlinear constraints
%  m                 = number of binding nonlinear constraints, plus one
%  tmax,LB,UB        = LP variables used for finding feasible descent direction
%  f,A,Aeq,b,beq     = LP variables used for finding feasible descent direction
%  OptOptions        = options used during LP solve
%  X,FX              = solution of LP solve
%  y                 = feasible descent direction (y = part of X)
%  exitflag          = error flag for bad LP problem
%  feasibleMeshPoint = flag indicating a feasible mesh point has been found
%  z                 = closest mesh point to pcenter + y
%  I                 = identity matrix
%  Ax                = product of Omega.A and iterate.x
%  feasible          = flag for indicating if poll center is feasible
%  lactive           = flags for numerically active lower bounds of Ax
%  uactive           = flags for numerically active upper bounds of Ax
%  active            = flags for numerically active bound/linear constraints
%  [B,N]             = tangent cone generating directions
%  Bgood,Bbad        = partitioning of B to deal with unavailable partials
%  Ngood,Nbad        = partitioning of N to deal with unavailable partials
%  gradcB            = least squares approx of NLP tangent cone generators
%  DN                = descent vectors in N
%  NN                = [N, -N]
%  descent           = descent vectors in I
%  dhat              = sum of the chosen ascent and decent directions
%*******************************************************************************
function D = gradDirections(iterate,Omega,strategy,Filter,delta,tol,scale)

% Determine type of pruning vector
switch strategy
    case 'Gradient_3n_L1'
        p = 1;
    case 'Gradient_3n_L2'
        p = 2;
    case 'Gradient_3n_LInf'
        p = Inf;
    case 'Gradient_3n2n'
        p = 2;
    otherwise
        error('mads:poll:choice','Invalid Poll strategy type (gradDirections).')
end

% Compute the descent directions that will be used for pruning
infGrad = ~isfinite(iterate.gradf);
g = iterate.gradf;
g(find(infGrad)) = 0;
d = -dVector(g.*scale,p);
descent = diag(-dVector(g,Inf));

% Determine the index of binding and violated constraints
B = find(iterate.c > -tol);
if (~isempty(B))
    
    %  If hmax is exceeded, then use d = -gradh
    if (iterate.h > Filter.hmax)
        infGrad = [~isfinite(iterate.gradh) infGrad];
        g = iterate.gradh;
        g(find(infGrad)) = 0;
        d = -dVector(g.*scale,p);
        descent = diag(-dVector(g,Inf));
        
        % If hmax is not exceeded, solve LP to get descent direction in f and cB
    else
        tmax = delta^2;
        m    = length(B)+1;
        f    = [-1; zeros(iterate.n,1); zeros(m,1)];
        A    = [ ones(m,1), zeros(m,iterate.n), -eye(m)];
        Aeq  = [zeros(m,1), [iterate.gradf'; iterate.gradc(:,B)'], eye(m)];
        b    =  zeros(m,1);
        beq  =  zeros(m,1);
        LB   = [0;   -Inf*ones(iterate.n,1); zeros(m,1)];
        UB   = [tmax; Inf*ones(m + iterate.n,1) ];
        OptOptions = optimset('Display','off');
        [X,FX,exitflag] = linprog(f,A,b,Aeq,beq,LB,UB,[],OptOptions);
        if (exitflag < 0), D = []; return, end
        y    = X(2:iterate.n+1);
        
        % Modify descent direction so that the resulting point is on the mesh
        m = 1;
        feasibleMeshPoint = 0;
        while (~feasibleMeshPoint)
            m = m + 1;
            z = iterate.x + delta*round((m*y - iterate.x)/delta);
            feasibleMeshPoint = iterate.gradf'*z<0 & all(iterate.gradc(:,B)'*z<0);
            if (m > 1e+3)
                error('mads:poll:badGradient', ...
                    'No feasible descent direction found (gradDirections).');
            end
        end
        d = z;
    end
end

% Add tangent cone generators and compensate for unavailable derivatives
I = eye(iterate.n);
if (isequal(Omega.A,I))
    Ax = Omega.A*iterate.x;
    lactive = activeConstraints( Ax, Omega.l, Omega.u,tol);
    uactive = activeConstraints(-Ax,-Omega.u,-Omega.l,tol);
    active  = lactive | uactive;
    Bgood   = find(~infGrad(:,end) & xor(lactive,uactive));
    Ngood   = find(~infGrad(:,end) & ~active);
    Bbad    = find( infGrad(:,end) &  active);
    Nbad    = find( infGrad(:,end) & ~active);
    B       = [descent(:,Bgood) I(:,Bbad) -I(:,Bbad)];
    DN      = I(:,Nbad);
    dhat    = -sum([-descent(:,Ngood) DN],2);
else
    [B,N]   = getTangentCone(iterate.x,Omega,tol,scale);
    Bgood   = find(B'*infGrad(:,end) == 0 & B'*g < 0);
    Ngood   = find(N'*infGrad(:,end) == 0);
    Bbad    = setdiff(1:size(B,2),Bgood);
    Nbad    = setdiff(1:size(N,2),Ngood);
    B       = [B(:,Bgood) B(:,Bbad) -B(:,Bbad)];
    DN      = N(:,Nbad);
    NN      = [N -N];
    dhat    = -sum([NN(:,NN'*g > 0) DN],2);
end
if (~any(dhat) || isequal(d(:,end),dhat)), dhat = []; end
if ~strcmp(strategy,'Gradient_3n2n'), descent = []; end
D = [scaleDirections(d,scale), dhat, scaleDirections(descent,scale),B,DN];
D(:,~any(D)) = [];
return

%*******************************************************************************
% getTangentCone:  Compute and return tangent cone generators.
% ------------------------------------------------------------------------------
% Acknowledgement:  This mathematics for this algorithm were developed by 
%    Robert Michael Lewis & Virginia Torczon (College of William & Mary)
% ------------------------------------------------------------------------------
% Called by: standDirections, gradDirections
% Calls:     activeConstraints, removeRedundancy
% VARIABLES:
%  B       = V*inv(V'*V)
%  N       = vectors than span the null space of V'
%  x       = current iterate continuous variables
%  Omega   = structure defining feasible region Omega = {x: l<=A*x<=u}
%    .A    =   matrix of linear constraint coefficients
%    .l    =   vector of constraint lower bounds
%    .u    =   vector of constraint upper bounds
%  tol     = tolerance within which constraints are considered active
%  scale   = scale factors for scaling the appropriate Poll directions
%  n       = number of continuous variables
%  I       = identity matrix
%  Ax      = A*x
%  V       = matrix whose columns are normal to active linear constraints
%  lactive = indices of linear constaints active at their lower bound
%  uactive = indices of linear constaints active at their upper bound
%  active  = indices of reformulated linear constraints
%  A       = alternative storage of linear constraint coefficients
%  b       = alternative storage of linear constraint bounds
%  [Q,R]   = QR factorization of V
%*******************************************************************************
function [B,N] = getTangentCone(x,Omega,tol,scale)

% Initialize variables
n = length(x);
I = eye(n);
Ax = Omega.A*x;
tolFinal = tol^2/2;

% Compute V
V = zeros(n,1);
while (rank(V) ~= min(size(V)))
    if (tol > tolFinal)
        lactive = activeConstraints( Ax,  Omega.l,  Omega.u, tol);
        uactive = activeConstraints(-Ax, -Omega.u, -Omega.l, tol);
        V = [ Omega.A(find(uactive),:)', -Omega.A(find(lactive),:)' ];
        tol = tol/2;
    else
        A = [Omega.A; -Omega.A]; 
        b = [Omega.u; -Omega.l];
        Ax = A*x;
        active = activeConstraints(Ax, [Omega.l; -Omega.u], b, tol);
        A = removeRedundancy(A, b, x, active);
        V = A';
        if rank(V) ~= min(size(V))
            error('mads:constraints:degenerate', ...
                ['Unable to remove degenerate nonredundant constraints ',...
                    '(getTangentCone).']);
        end
    end
end

% Compute B and N
if (isempty(V))
    B = zeros(n,1);
    B(:,1:end) = [];
    N = scaleDirections(I,scale);
else
    [Q,R] = qr(V,0);
    B = Q/R';
    N = I - B*V';
end

% Zero out elements that are very small already
N(abs(N) < 10*eps) = 0;
B(abs(B) < 10*eps) = 0;
return

%*******************************************************************************
% scaleDirections:  Scale each column by a vector of scale factors.
% ------------------------------------------------------------------------------
% Called by: standDirections, madsDirections, gradDirections
% VARIABLES:
%  Y     = output matrix
%  X     = input matrix
%  scale = vector of scale factors
%*******************************************************************************
function Y = scaleDirections(X,scale)
if (isempty(X) || all(scale == 1))
    Y = X;
else
    Y = diag(scale)*X;
end
return

%*******************************************************************************
% activeConstraints:  Check which constraint bounds are numerically active.
% ------------------------------------------------------------------------------
% Called by: standDirections, gradDirections, getTangentCone
% VARIABLES:
%  active = logicals indicating which constraints is numerically active
%  x      = point in R^n to be tested 
%  [a,b]  = upper and lower bounds on x
%  tol    = error tolerance for e-active constraints
%  infa   = indices of x having no lower bound
%  infb   = indices of x having no upper bound, but having a lower bound
%  finite = indices of x having finite lower and uppr bounds
%*******************************************************************************
function active = activeConstraints(x,a,b,tol)

% Set status of each continuous variable
infa   = find(isinf(a));
infb   = find(isinf(b) & ~isinf(a));
finite = isfinite(a.*b);
one    = ones(length(infb),min(1,length(infb)));

% Determine status of each variable
active         = zeros(length(x),1);
scale          = max([abs(x(infb)),abs(a(infb)),one],[],2);
active(infb)   = (x(infb)  -a(infb)   <= tol*scale);
scale          = max(b(finite)-a(finite),1);
active(finite) = (x(finite)-a(finite) <= tol*scale);
return

%*******************************************************************************
% removeRedundancy:  Removes rows of Ax <= b that are redundant.
% ------------------------------------------------------------------------------
% This algorithm is courtesy of Olga A. Brezhneva (Miami University, Ohio)
% ------------------------------------------------------------------------------
% Called by: getTangentCone
% VARIABLES:
%  AA     = matrix A with only nonredundant e-active constraints left
%  bb     = vector b with only nonredundant e-active constraints left
%  A      = coefficient matrix
%  b      = right-hand side
%  x      = current iterate continuous variables
%  active = indices of the active constraints (rows) in A
%  ind    = vector of indices
%  bb     = b vector with only nonredundant e-active constraints left
%  m      = length of bb
%  res    = residual of current constraint
%  px     = projection of x onto the current constraint boundary
%  pres   = vector of residuals with respect to px
%  NR     = logical indicating non-redundancy status
%  xval   = optimal x with respect to the LP solver
%  fval   = LP objective function value at xval
%  i      = constraint counter
%*******************************************************************************
function [AA,bb] = removeRedundancy(A,b,x,active)

% Form e-active constraints
ind = find(active);
AA  = A(ind,:);
bb  = b(ind);

% Begin loop through constraints
i = 1;
m = length(bb);
while i <= m
    
    % Strategy 1: verification of nonredundancy
    res  = bb(i) - AA(i,:)*x;
    px   = x + AA(i,:)'*res;
    pres = bb - AA*px;
    pres(i) = 1;
    NR = all(pres > 0);
    
    % Strategy 2: solve LP problem
    if ~NR
        ind    = 1:length(bb);
        ind(i) = [];
        warning('off','all');
        [xval,fval] = linprog(-AA(i,:)',AA(ind,:),bb(ind),[],[],[],[],[], ...
            optimset('Display','off'));
        warning('on','all');
        NR = (-fval > bb(i));
    end
    
    % Delete constraint if redundant
    if NR
        i = i + 1;
    else
        AA(i,:) = [];
        bb(i)   = [];
        m = length(bb);
    end
end
return

%*******************************************************************************
% dVector:  Returns the appropriate 3^n vector approximating v.
% ------------------------------------------------------------------------------
% Called by: gradDirections
% VARIABLES:
%  w = L-1, closest L-2, or L-Inf {-1,0,1}^n approximation to v
%  v = input vector
%  p = norm of the gradient
%*******************************************************************************
function w = dVector(v,p)
switch p
    case {1}
        w = (abs(v) == max(abs(v))).*sign(v);
    case {2}
        w = v/norm(v,inf);
        w(~isfinite(w) | abs(w) < tan(pi/8)) = 0;
        w = sign(w);
    case {inf}
        w = sign(v);
end
return

%*******************************************************************************
% getPollOrder:  Returns the order in which to poll from the poll set.
% ------------------------------------------------------------------------------
% Called by: poll
% VARIABLES:
%  order     = resulting Poll order
%  orderCode = string indicating user choice
%  D         = set of Poll directions
%  P         = the Poll set
%  RunData   = structure of Run statistics and information
%    .goodD  =   previous direction of the last successful iteration
%    .porder =   poll order from the previous successful iteration
%  n         = dimension of the optimization problem
%  nD        = number of Poll directions
%  norms     = norms of the columns of D
%  angles    = angles between .goodD and the columns of D
%  sf,sc,sh  = surrogate f-, c-, and h-values
%*******************************************************************************
function order = getPollOrder(orderCode,D,P,RunData,surrogate)

% get the dimensions of D
[n,nD] = size(D);

% Construct the order according to user choice
switch orderCode
    case 'Consecutive'
        order = 1:nD;
    case 'Alternating'
        switch nD
            case {n+1}
                order = [nD, 1:nD-1];
            case {2*n}
                order = reshape([1:n; n+1:2*n],1,nD);
            otherwise
                order = 1:nD;
        end
    case 'Random'
        order = randperm(nD);
    case {'Dynamic','DynamicRanked'}
        order = 1:nD;
        if isfield(RunData,'porder')
            norms = zeros(nD,1);
            for k = 1:nD
                norms(k) = norm(D(:,k));
            end
            angles = (D'*RunData.goodD)./norms;
            [y,order] = sort(-angles);
            if strcmp(orderCode,'Dynamic')
                k = find(RunData.porder == order(1));
                order = [order(1) RunData.porder(1:k-1) RunData.porder(k+1:end)];
            end
        end
    case 'Surrogate'
        sf = feval(surrogate.evaluator,[P.x]',surrogate.f);
        if isfield(surrogate,'c')
            for k = 1:length(surrogate.c)
                sc(:,k) = feval(surrogate.evaluator,[P.x]',surrogate.c(k));
            end
            sh = sum((sc > 0).*sc.*sc,2);
        end
        if (RunData.pcenter.h > 0)
            [y,order] = sort(sh);
        else
            [y,order] = sort(sf);
        end
    case 'Custom'
        order = RunData.porder;
        if (length(order) ~= nD)
            error('mads:poll:dim', ...
                ['Poll directions and Poll orders have incompatible dimensions ', ...
                    '(getPollOrder).']);
        end
end
return

%*******************************************************************************
% FUNCTIONS CALLED BY MADS UPDATE ROUTINE
%*******************************************************************************
%*******************************************************************************
% getPollCenter: Gets the POLL center based on CenterCode.
% ------------------------------------------------------------------------------
% Called by: processInput, processOutput, search, mvpPoll, update
% VARIABLES:
%  pcenter      = chosen POLL center
%  iCache       = Cache index of poll center
%  centerCode   = code indicating which point to POLL around
%                 (0 = current solution; otherwise = Filter(centerCode))
%  Filter       = sorted indices of Cache iterates in the filter
%  Cache        = collection of all previous iterates
%  closest      = min(centerCode, filter size)
%  noFilter     = flag indicating no infeasible iterates have been found
%  noFeasible   = flag indicating no feasible iterates have been found
%  filterPoints = indices of valid filter points
%*******************************************************************************
function [pcenter,iCache] = getPollCenter(centerCode,Filter,Cache)

% Determine status of filter vectors.  Error if both are empty.
noFilter   = isempty(Filter.F);
noFeasible = isempty(Filter.feasible);
if (noFeasible && noFilter)
    error('mads:pollCenter:badFilter', ...
        'No suitable iterate to poll around (PollCenter).');
end

% Want current best feasible point.  If none, pick least infeasible point.
if (centerCode == 0)
    if (noFeasible)
        iCache = Filter.F(1);
    else
        iCache = Filter.feasible(1);
    end
    
    % Want particular filter point.  If none, pick current best feasible point.
else
    if (noFilter)
        iCache = Filter.feasible(1);
    else
        filterPoints = find([Cache.iterate(Filter.F).h] < Filter.hmax);
        if ~isempty(filterPoints)
            closest = min(centerCode,filterPoints(end));
        else
            closest = 1;
        end
        iCache = Filter.F(closest);
    end
end

pcenter = Cache.iterate(iCache);
return;

%*******************************************************************************
% updateOmega:  Update the set of linear constraints .l <= .A*x <= .u.
% ------------------------------------------------------------------------------
% Parts of this routine are courtesy of Dennis & Schnabel textbook
% ------------------------------------------------------------------------------
% Called by: processInput, search, mvpPoll, update, evalPointSet
% Calls:     getpID, < Omega File >
% VARIABLES:
%  Omega               = structure describing (linear) feasible region
%    .A                =   coefficient matrix
%    .l                =   vector of lower bounds
%    .u                =   vector of upper bounds
%    .plist            =   list of allowable categorical variable values
%  newScale            = scale factors for constructing mesh
%  Problem             = structure of problem data
%    .File.O           =   name of file defining feasible region Omega
%    .isMVP            =   flag indicating problem is an MVP problem
%  Options             = structure of MADS options
%    .removeRedundancy =   turns on/off removal of redundant constraints
%    .scale            =   base for logarithmic scaling
%  iterate             = current iterate
%    .n                =   dimension of the continuous variable space
%    .x                =   vector of continuous variable values
%    .p                =   cell array of categorical variables
%  n                   = number of continuous variables of iterate
%  m                   = number of linear constraints
%  [a,b]               = temporary storage of variable bounds
%  scaleA              = vector of scales for scaling linear constraints
%  AA,bb               = reformulation of linear constraints
%  Redundant           = flags indicating redundant linear constraints
%  ind                 = indices of nonredundant linear constraints
%  indcopy             = temporary copy of ind
%  xval,fval           = x-value and f-value of LP solution
%**************************************************************************
function [Omega,newScale] = updateOmega(Problem,Options,iterate)

% Get parameters that define the linear constraints
n = length(iterate.x);
if exist(Problem.File.O,'file')
    switch nargin(Problem.File.O) + 2*Problem.isMVP
        case {1}
            [A,l,u] = feval(Problem.File.O,n);          
        case {2}
            error('mads:constraints:badFileFormat', ...
                'Non-MVP Omega file cannot have a second argument (updateOmega).');
        case {3}
            [A,l,u,plist] = feval(Problem.File.O,n);    
        case {4}
            [A,l,u,plist] = feval(Problem.File.O,n,iterate.p);    
    end
else
    if (Problem.isMVP)
        error('mads:constraints:file', ...
            ['MVP problem missing its Omega file, ',Problem.File.O, ...
                ' (updateOmega).']);
    else
        [A,l,u] = deal(eye(n),-Inf*ones(n,1),Inf*ones(n,1));
    end
end
[m,n] = size(A);

% Error check linear constraint parameters
errorCheck = ~isempty(find(l - u > 0)) || size(l,1) ~= m || ...
    all(size(l) ~= size(u))    || size(u,1) ~= m;
if (errorCheck),
    error('mads:constraints:dim','Invalid constraint dimensions (updateOmega).');
end

% Scale the linear constraints by their inf-norms and form Omega structure
a = l(1:n);
b = u(1:n);
scaleA = max(abs(A),[],2);
A = diag(1./scaleA)*A;
l = l./scaleA;
u = u./scaleA;

%% Detect redundant linear constraints (due in part to Olga Brezhneva)
ind = 1:m;
if (Options.removeRedundancy)
    AA = [-A; A];
    bb = [-l; u;];
    Redundant = ~isfinite(bb)';
    ind = find(isfinite(bb));
    for i = 1:length(ind)
        indcopy = ind;
        indcopy(i) = [];
        warning('off','all');
        [xval,fval] = linprog(-AA(ind(i),:)',AA(indcopy,:),bb(indcopy), ...
            [],[],[],[],[],optimset('Display','off'));
        warning('on','all');
        if ~isempty(fval)
            Redundant(ind(i)) = (-fval <= bb(ind(i)));
        else
            error('mads:constraints:LPerror', ...
                ['linprog solver failed. Recheck LP formulation for removing ',...
                    'redundant linear constraints (updateOmega).']);
        end
    end
    
    % Delete redundant linear constraints
    ind = all(reshape(Redundant,m,2)');
    A(ind,:)   = [];
    l(ind)     = [];
    u(ind)     = [];
end

% Form Omega structure and get new scale factors
Omega = struct('A',full(A),'l',l,'u',u);
if (Problem.isMVP)
    Omega.plist = plist;
    newScale = getScaleFactors(Options.scale,iterate,a,b,Problem.nameCache, ...
        plist);
else
    newScale = getScaleFactors(Options.scale,iterate,a,b,Problem.nameCache);
end
return

%*******************************************************************************
% getScaleFactors:  Get scale factors for scaling directions.
% ------------------------------------------------------------------------------
% Called by: UpdateOmega
% VARIABLES:
%  newScale  = vector of scale factors for scaling of directions
%  scale     = base of logarithmic scaling (0 = no scale)
%  iterate   = current poll center
%  [a,b]     = vectors of lower and upper bounds on continuous variables
%  isMVP     = flag indicating problem is an MVP problem
%  nameCache = name of the base workspace Cache variable
%  typx      = "typical" values used for scaling if bounds are insufficient
%  pID       = ID for the categorical variable values
%  ind       = index used to match Cache points with the same pID
%  zeroInf   = flags for zero or infinite values in [a,b]
%  bad       = indices of variables with a and b both zero or infinite
%  good      = indices of variables with a and b both finite and nonzero
%  lhalf     = indices of variables with finite nonzero a and zero or infinite b
%  rhalf     = indices of variables with finite nonzero b and zero or infinite a
%  expScale  = exponents of the approximate ranges of the variables
%  baseScale = factor for coverting from base b to base 10 logarithm
%*******************************************************************************
function newScale = getScaleFactors(scale,iterate,a,b,nameCache,varargin)

% Compute scale factors, if selected
n = length(iterate.x);
if (scale)
    
    % Compute "typical" value of x for scaling
    isMVP = (nargin > 5);
    Cache = getappdata(0,nameCache);
    if isempty(Cache.iterate)
        typx = iterate.x;
    else
        if (isMVP)
            pID = getpID(iterate.p,varargin{:});
            ind = find(Cache.pID(1:Cache.size) == pID);
            if isempty(ind)
                typx = iterate.x;
            else
                typx = Cache.iterate(ind(1)).x;
            end
        else
            typx = Cache.iterate(1).x;
        end
    end
    
    % Flag zero or infinite bounds, and categorize each variable's bounds
    zeroInf = [ (~isfinite(a)|abs(a)<=eps), (~isfinite(b)|abs(b)<=0)];
    good    = find(~any(zeroInf,2));
    lhalf   = find(~zeroInf(:,1)  &  zeroInf(:,2));
    rhalf   = find( zeroInf(:,1)  & ~zeroInf(:,2));
    bad     = find(all(zeroInf,2) &  typx);
    verybad = find(all(zeroInf,2) & ~typx);
    
    % Scale variables with finite nonzero bounds IAW D&S, ch 7 and pp. 278-9
    baseScale       = 1/log10(scale);
    expScale        = ones(n,1);
    expScale(good)  = baseScale * (log10(abs(a(good)))+log10(abs(b(good))))/2;
    expScale(lhalf) = baseScale *  log10(abs(a(lhalf)));
    expScale(rhalf) = baseScale *  log10(abs(b(rhalf)));
    expScale(bad)   = baseScale *  log10(abs(typx(bad)));
    expScale(verybad) = 0;
    newScale = scale.^round(expScale);
else
    newScale = ones(n,1);
end
return

%*******************************************************************************
% FUNCTIONS FOR EVALUATING FUNCTIONS AND UPDATING CACHE DATA
%*******************************************************************************
%*******************************************************************************
% evalPointSet:  Evaluates the objective and constraint functions.
% ------------------------------------------------------------------------------
% Called by: processInput, search, poll, mvpPoll
% Calls:     evalRSPointSet, updateOmega, inOmega, getpID, isCacheHit, 
%            evalFunction, updateCache, updateFilter
% VARIABLES:
%  P              = set of iterates to be evaluated
%  success        = logical indicating a successful iteration
%  Filter         = updated filter
%  RS             = structure of R&S parameters (unused if not stochastic)
%  Problem        = structure containing problem parameters
%    .File.O      =   name of Omega file
%    .fType       =   type of functions file (C=C, F=FORTRAN, M=MATLAB)
%    .isMVP       =   logical indicating if problem is n MVP
%    .Omega       =   structure defining feasible region {x: l<=Ax<=u}
%      .A         =     coefficient matrix
%      .l         =     vector of lower bounds
%      .u         =     vector of upper bounds
%      .plist     =     list of valid categorical variables
%    .nameCache   =   name of the base workspace Cache variable
%  full           = flag indicating if all points in P will be evaluated
%  Options        = structure of user options
%    .computeGrad =   flag for computing gradients
%  nP             = number of iterates to be evaluated
%  pID            = categorical variable value ID
%  xnorm          = 1-norm of iterate being evaluated
%  hit            = flag indicating a Cache hit
%  unfiltered     = logical indicating an iterate is unfiltered
%*******************************************************************************
function [P,success,Filter,RS] = evalPointSet(Problem,P,full,Options,Filter,...
    RS,add2sur)

success = 0;
if (isempty(P)), return; end

% Initialize output fields
nP = length(P);
for k = 1:nP
    P(k).n = length(P(k).x);
end
[P.f,P.c,P.h] = deal([]);

% Call R&S version if R&S problem
if Options.runStochastic
    [P,success,Filter,RS] = evalRSPointSet(Problem,P,Options,Filter,RS,add2sur);
    return
end

if (Options.computeGrad)
    [P.gradf,P.gradc,P.gradh] = deal([]);
end

% Evaluate each iterate, as appropriate
for k = 1:nP
    
    % update MVP Omega parameters, as is necessary
    if Problem.adaptOmega
        [Problem.Omega] = updateOmega(Problem,Options,P(k));
    end
    
    % Test if iterate lies in Omega
    if (inOmega(P(k),Problem.Omega))
        
        % Give trial point unique IDs to expedite Cache searching
        pID = 0;
        if Problem.isMVP
            pID = getpID(P(k).p,Problem.Omega.plist);
        end
        xnorm = norm(P(k).x,inf);
        
        % Test to see if the iterate is in the Cache
        Cache = getappdata(0,Problem.nameCache);
        hit = isCacheHit(Cache,P(k),xnorm,pID,50000);
        Cache.nHits = Cache.nHits + hit;
        setappdata(0,Problem.nameCache,Cache);
        
        % Evaluate function, store in Cache, and update Cache and filter
        if (~hit)
            P(k) = evalFunction(Problem.File.F,Problem.fType,P(k), ...
                Options.computeGrad);
            Cache = updateCache(Cache,Problem,P(k),xnorm,pID,1);
            setappdata(0,Problem.nameCache,Cache);
            if isfinite(P(k).h)
                [unfiltered,Filter] = updateFilter(P(k),Filter,Cache);
            else
                unfiltered = 0;
            end
            if (unfiltered), success = k; end
            if (add2sur == 2) || (add2sur == 1 && unfiltered)
                Cache = getappdata(0,Problem.nameCache);
                Cache.isSurrPoint(Cache.size) = 1;
                setappdata(0,Problem.nameCache,Cache);
            end
            if (success && ~full), break; end
        end
    end
end
return

%*******************************************************************************
% getpID:  Get unique ID based on values of categorical variables.
% ------------------------------------------------------------------------------
% Called by: evalPointSet
% VARIABLES:
%  pID   = unique ID
%  p     = cell array of categorical variables
%  plist = cell array of allowable values for each categorical variable
%  m     = vector of plist lengths for each categorical variable
%  ind   = index of plist that matches the current categorical variable 
%*******************************************************************************
function pID = getpID(p,plist)

np  = length(p);
m   = zeros(np,1);
pID = zeros(np,1);
for k = 1:np
    m(k) = length(plist{k});
    if ischar(p{k})
        ind = find(strcmp(p{k},plist{k}));
    else
        ind = find([plist{k}{:}] == p{k});
    end
    pID(k) = ind(1);
end
pID = pID'*[1; cumprod(m(1:end-1)+1)];
return

%*******************************************************************************
% inOmega:  Test to see if point x satisfies l <= A*x <=u.
% ------------------------------------------------------------------------------
% Called by: makeFeasible, evalPointSet
% VARIABLES:
%  pass     = logical indicating is x satisfies l <= A*x <=u
%  iterate  = current iterate
%    .x     =   vector of continuous variables
%    .p     =   cell array of categorical variables
%  Omega    = feasible space with respect to bound and linear constraints
%    .A     = coefficient matrix for linear and bound constraints
%    .l     = lower bound vector for linear and bound constraints
%    .u     = upper bound vector for linear and bound constraints
%    .plist = cell array of lists of allowed categorical values in p
%  Ax       = .A*x
%*******************************************************************************
function pass = inOmega(iterate,Omega)

% Check for errors in categorical variables or Omega construction
for k = 1:length(iterate.p)
    if (ischar(iterate.p{k}))
        pass = any(strcmp(iterate.p{k},Omega.plist{k}));
    else
        pass = ~isempty(find([Omega.plist{k}{:}] == iterate.p{k}));
    end
    if (~pass), return; end
end

% Check if iterate is in Omega
Ax = Omega.A*iterate.x;
pass = isempty(find(Omega.l > Ax)) && isempty(find(Omega.u < Ax));
return;

%*******************************************************************************
% isCacheHit:  Test to see if point has already been evaluated previously.
% ------------------------------------------------------------------------------
% Called by: evalPointSet
% VARIABLES:
%  hit        = logical indicating that the iterate matches a Cache point
%  Cache      = array of previously computed iterates
%    .iterate =   vector of iterates
%    .tol     =   tolerance for declaring a Cache hit w.r.t iterate.x
%    .xnorm   =   vector of 1-norms of iterates in the Cache
%    .pID     =   vector of IDs od categorical variable values
%  iterate    = current point
%    .x       =   coordinates of continuous  variables
%    .p       =   coordinates of categorical variables
%  xnorm      = 1-norm of iterate.x
%  pID        = ID of categorical variable values
%  looksize   = how far back into the Cache will be searched
%  ind        = indices of Cache points that may match the current point
%  xvalues    = temporary storage of Cache continuous variable values
%*******************************************************************************
function hit = isCacheHit(Cache,iterate,xnorm,pID,looksize)

% Initialize Cache search
hit = 0;
ind = max(1,Cache.size-looksize):Cache.size;

% Screen out iterates without sufficiently close x-norms (see Apostol, p49)
ind = ind((abs(Cache.xnorm(ind) - xnorm) < Cache.tol));
if isempty(ind), return, end

% Screen out iterates with different categorical variable values
if pID ~= 0
    ind = ind(Cache.pID(ind) == pID);
    if isempty(ind), return, end
end

% Finish search
xvalues = [Cache.iterate(ind).x];
ind = ind(max(abs(xvalues - repmat(iterate.x,1,length(ind))),[],1) < Cache.tol);
hit = ~isempty(ind);
return

%*******************************************************************************
% evalFunction:  Evaluates objective and constraint functions at a point.
% ------------------------------------------------------------------------------
% Called by: evalPointSet
% Calls: < User F >
% VARIABLES:
%  f        = name of optimization problem functions file
%  ftype    = type of functions file (F=FORTRAN, C=C/C++, M=MATLAB)
%  iterate  = current iterate
%    .x     =   coordinates of iterate
%    .p     =   values of categorical variables
%    .f     =   objective function value of iterate
%    .c     =   constraint function values of iterate
%    .h     =   constraint violation function value of iterate
%    .gradf =   f'(x) at current iterate x
%    .gradc =   c'(x) at current iterate x
%    .gradh =   h'(x) at current iterate x
%  grad     = flag for computing gradients, if available
%  pchar    = categorical variables stored as character strings
%  pint     = categorical variables stored as integers
%  cflag    = flag for sorting out char and int categorical variables
%  nc       = number of nonlinear constraints
%  badValue = the value used when f(x) is complex (Inf)
%  cplus    = vector of nonlinear constraint violations
%*******************************************************************************
function iterate = evalFunction(f,ftype,iterate,grad)

switch upper(ftype)
    
    % Process a compiled Fortran functions file   
    case {'F','C'}
        if (isempty(iterate.p))
            pchar = 0;
            pint  = 0;
        else
            cflag = zeros(length(iterate.p),1);
            for k = 1:length(iterate.p)
                cflag(k) = ischar(iterate.p{k});
            end
            pchar = [iterate.p{find( cflag)}];
            pint  = [iterate.p{find(~cflag)}];
        end
        
        [nc,iterate.f,c,gradf,gradc] = feval(f,iterate.x,pint,pchar);
        iterate.c = c(1:nc);
        if (grad)
            iterate.gradf = gradf;
            iterate.gradc = gradc(:,1:nc);
        end
        
        % Process a Matlab functions file   
    case {'M'}
        switch nargin(f) + 2*grad;
            case {1}
                if (nargout(f) == 1)
                    iterate.f = feval(f,iterate.x);
                    iterate.c = [];
                else
                    [iterate.f,iterate.c] = feval(f,iterate.x);
                end
            case {2}
                if (nargout(f) == 1)
                    iterate.f = feval(f,iterate.x,iterate.p);
                    iterate.c = [];
                else
                    [iterate.f,iterate.c] = feval(f,iterate.x,iterate.p);
                end
            case {3}
                [iterate.f,iterate.c,iterate.gradf,iterate.gradc] = feval(f,iterate.x);
            case {4}
                [iterate.f,iterate.c,iterate.gradf,iterate.gradc] = feval(f,iterate.x, ...
                    iterate.p);
            otherwise
                error('mads:function:badVariables', ...
                    'Error in FunctionEval logical input variables (functionEval).');
        end
    otherwise
        error('mads:function:badType','Invalid Function Type (functionEval).');
end

% Prevent complex-valued functions
badValue = realmax;
if (~isreal(iterate.f)), iterate.f = badValue; end
iterate.c(~isreal(iterate.c)) = badValue;
if (grad)
    iterate.gradf(~isreal(iterate.gradf)) = badValue;
    iterate.gradc(~isreal(iterate.gradc)) = badValue;
end

% Compute h(x) (and h'(x) if available)
cplus = (iterate.c > 0).*iterate.c;
iterate.h = norm(cplus)^2;
if (grad)
    iterate.gradh = 2*iterate.gradc*cplus';
end

return;

%*******************************************************************************
% updateCache:  Update the Cache with the current iterate.
% ------------------------------------------------------------------------------
% Called by: EvalPointSet
% VARIABLES:
%  Cache      = structure containing data on all past iterates
%    .iterate =   the previously computed iterates
%    .size    =   current number of iterates in the Cache
%    .tol     =   tolerance for determining if iterate was in Cache
%    .xnorm   =   vector of 1-norms of previously computed iterates
%    .pID     =   vector of categorical variable ID numbers
%    .bfp     =   vector of best feasible points
%    .lip     =   vector of least infeasible points
%  Problem    = structure containing problem file names and parameters
%    .isMVP   =   logical indicating if problem is an MVP
%    .maxNc   =   maximum number of constraints to be added
%    .maxNx   =   maximum number of continuous variables of any iterates
%    .maxNp   =   maximum number of categorical variables of any iterates
%  iterate    = current iterate
%    .x       =   coordinates of iterate
%    .p       =   values of categorical variables
%    .n       =   number of continuous variables
%    .f       =   objective function value of iterate
%    .c       =   constraint function values of iterate
%    .h       =   constraint violation function value of iterate
%    .gradf   =   f'(x) at current iterate x
%    .gradc   =   c'(x) at current iterate x
%    .gradh   =   h'(x) at current iterate x
%  xnorm      = 1-norm of the curent iterate
%  pID        = categorical variable value ID of the current iterate
%  nFunc      = number of function evaluations at point (=1, if not stochastic)
%  n          = number of iterates for which memory is to be allocated
%  ind        = index variable for additional chunks of memory
%  p          = temporary storage
%*******************************************************************************
function Cache = updateCache(Cache,Problem,iterate,xnorm,pID,nFunc)

% Allocate extra chunks of memory if more is required 
if (Cache.size >= max(1024,length(Cache.iterate)))
    n = 1024;
    ind = length(Cache.iterate) + (1:n);
    if (Problem.isMVP)
        [p{1:Problem.maxNp}] = deal('            ');
    else
        p = {};   
    end
    [Cache.iterate(ind).x] = deal(zeros(Problem.maxNx,1));
    [Cache.iterate(ind).p] = deal(p);
    [Cache.iterate(ind).n] = deal(0);
    [Cache.iterate(ind).f] = deal(0);
    [Cache.iterate(ind).c] = deal(zeros(Problem.maxNc,1));
    [Cache.iterate(ind).h] = deal(0);
    Cache.xnorm(ind)       = deal(0);
    Cache.pID(ind)         = deal(0);
    Cache.bfp(ind)         = deal(0);
    Cache.lip(:,ind)       = deal(0);
    Cache.nFunc(ind)       = deal(0);
    Cache.isSurrPoint(ind) = deal(0);
    if (isfield(iterate,'gradf'))
        [Cache.iterate(ind).gradf] = deal(zeros(Problem.maxNx,1));
        [Cache.iterate(ind).gradc] = deal(zeros(Problem.maxNc,Problem.maxNx));
        [Cache.iterate(ind).gradh] = deal(zeros(Problem.maxNx,1));
    end
end

% Update Cache (.isSurrPoint is set to zero and updated later)
Cache.size = Cache.size + 1;
Cache.iterate(Cache.size)     = iterate;
Cache.xnorm(Cache.size)       = xnorm;
Cache.pID(Cache.size)         = pID;
Cache.nFunc(Cache.size)       = nFunc;
Cache.isSurrPoint(Cache.size) = 0;

% Update the BFP and LIP function values in the Cache (used in plotting)
feasible = (iterate.h < Cache.Filter.hmin);
if (feasible)
    if (Cache.size == 1)
        Cache.bfp(Cache.size) = iterate.f;
        Cache.lip(:,Cache.size) = [Inf; Inf];
    else
        Cache.bfp(Cache.size) = min(Cache.bfp(Cache.size-1),iterate.f);
        Cache.lip(:,Cache.size) = Cache.lip(:,Cache.size-1);
    end
else
    if (Cache.size == 1)
        Cache.bfp(Cache.size) = Inf;
        Cache.lip(:,Cache.size) = [iterate.f; iterate.h];
    else
        Cache.bfp(Cache.size) = Cache.bfp(Cache.size-1);
        if (iterate.h < Cache.lip(Cache.size-1))
            Cache.lip(:,Cache.size) = [iterate.f; iterate.h];
        else
            Cache.lip(:,Cache.size) = Cache.lip(:,Cache.size-1);
        end
    end
end

return

%*******************************************************************************
% FUNCTIONS FOR UPDATING AND PLOTTING THE FILTER
%*******************************************************************************
%*******************************************************************************
% updateFilter:  Update the filter and solutions vectors, given in iterate.
% ------------------------------------------------------------------------------
% Called by: evalPointSet, mvpPoll
% Calls:     dominates, plotFilter
% VARIABLES:
%  unfiltered  = logical indicating that the iterate is unfiltered
%  Filter      = structure containing filter
%    .hmin     =   minimum allowed constraint violation of filter point
%    .hmax     =   maximum allowed constraint violation of filter point
%    .strict   =   flag indicating that hmax is strictly enforced
%    .plot     =   flag for displaying a real-time plot of the filter 
%    .feasible =   sorted indices of feasible Cache iterates
%    .F        =   sorted indices of Cache iterates in the filter
%  iterate     = current iterate
%    .x        =   coordinates of iterate
%    .p        =   values of categorical variables
%    .n        =   dimension of the continuous variables
%    .f        =   objective function value of iterate
%    .c        =   constraint function values of iterate
%    .h        =   constraint violation function value of iterate
%    .gradf    =   f'(x) at current iterate x
%    .gradc    =   c'(x) at current iterate x
%    .gradh    =   h'(x) at current iterate x
%  Cache       = collection of all previously evaluated iterates
%    .iterate  =   vector of iterates
%    .size     =   number of iterates
%  maxSolPts   = maximum number of feasible points to keep in filter
%  ind         = indices used in sorting
%  dominated   = flags indicating which filter points iterate dominates
%  infeasible  = flags indicating which filter points exceed hmax
%*******************************************************************************
function [unfiltered,Filter] = updateFilter(iterate,Filter,Cache)

maxSolPts = 10;
feasible  = iterate.h <= Filter.hmin;

% FEASIBLE CASE: Update sorted list of 10 best feasible solutions
if (feasible)
    fvalues = [Cache.iterate(Filter.feasible).f];
    ind = find(fvalues > iterate.f);
    if (isempty(ind))
        if (isempty(Filter.feasible)), Filter.feasible(1) = Cache.size; end
    else
        Filter.feasible = [Filter.feasible(1:ind(1)-1); Cache.size; ...
                Filter.feasible(ind(1):end)];
        Filter.feasible(maxSolPts:end) = [];
    end
    unfiltered = (Filter.feasible(1) == Cache.size);
    return;
end

% INFEASIBLE CASE

% Initialize hmax and test to see if point is unfiltered
if (~Filter.strict && isempty(Filter.F))
    hmax = 10*iterate.h;
elseif (~Filter.strict && Cache.iterate(Filter.F(1)).h >= Filter.hmax)
    hmax = 10*Cache.iterate(Filter.F(1)).h;
else
    hmax = Filter.hmax;
end
unfiltered = iterate.h<hmax && ~any(dominates(Cache.iterate(Filter.F),iterate));

% Delete newly dominated points, add the unfiltered point, and replot
if (unfiltered)
    dominated  = dominates(iterate,Cache.iterate(Filter.F));
    infeasible = [Cache.iterate(Filter.F).h] >= hmax;
    Filter.F(dominated | infeasible) = [];
    hvalues = [Cache.iterate(Filter.F).h];
    ind     = find(hvalues > iterate.h);
    if (isempty(ind))
        Filter.F(end+1) = Cache.size;
    else
        Filter.F = [Filter.F(1:ind(1)-1), Cache.size, Filter.F(ind(1):end)];
    end
    if (Filter.plot), plotFilter(Filter,hmax,Cache.iterate,maxSolPts); end
end
return;

%*******************************************************************************
% dominates:  Determines if x dominates y with respect to .f and .h values.
% ------------------------------------------------------------------------------
% Called by: updateFilter
% VARIABLES:
%  d    = logicals indicating whether or not x dominates y
%  x,y  = two iterates to be compared
%    .f =   objective function value
%    .h =   constraint violation function value
%*******************************************************************************
function d = dominates(x,y)
d = ~([y.f] < [x.f] | [y.h] < [x.h]);
return

%*******************************************************************************
% plotFilter:  Plot the filter.
% ------------------------------------------------------------------------------
% Called by: updateFilter, processInput
% VARIABLES:
%  Filter        = structure containing the filter
%    .F          = sorted indices of Cache iterates in the filter
%    .feasible   = sorted indices of feasible Cache iterates
%    .hmax       = maximum allowed constraint violation of filter point
%    .plothandle = plot handle
%  hmax          = the filter hmax for the plot only
%  Cache         = collection of all previously evaluated iterates
%    .f          =   objective function value
%    .h          =   constraint violation function value
%  nPoints       = number of filter points to plot
%  nFilter       = number of iterates in the filter
%  nSolutions    = number of feasible solutions
%  nPlot         = number of feasible solutions to plot
%  feasible_f    = vector of f-values of feasible points to be plotted
%  h             = vector of h-values of filter points to be plotted
%  f             = vector of f-values of filter points to be plotted
%  pctAxis       = percentage of extra space between axes and data
%  ax            = plot axis parameters
%*******************************************************************************
function plotFilter(Filter,hmax,Cache,nPoints)

if (isempty(Filter.F) || isempty(Filter.feasible)), return; end
nSolutions = length(Filter.feasible);
nPlot      = min(nPoints,nSolutions);

% Create feasible solutions vector and filter vectors for plotting
feasible_f = [Cache(Filter.feasible(:)).f];
[h,f] = stairs([Cache(Filter.F).h],[Cache(Filter.F).f]);
h = [h(1); h];
f = [max(Cache(Filter.feasible(nPlot)).f,Cache(Filter.F(1)).f); f];
h(end+1) = hmax;
f(end+1) = f(end);

% Create filter plot
axes(Filter.plothandle);
set([Filter.plothandle; get(Filter.plothandle,'Children')],'Visible','on');
plot(h,f,'k-',zeros(nPlot,1),feasible_f(1:nPlot),'b*');
pctAxis = 0.02;
ax = axis;
ax(1:2) = ax(1:2) + pctAxis*(ax(2)-ax(1))*[-1,1];
ax(3:4) = ax(3:4) + pctAxis*(ax(4)-ax(3))*[-1,1];
line([h(end),h(end)],[f(end) ax(3)]);
axis(ax);
title('Filter','FontWeight', 'bold', 'FontSize',11);
xlabel('Constraint Violation, h(x)', 'FontSize',12);
ylabel('Objective Function, f(x)  ', 'FontSize',12);
drawnow;

return;

%*******************************************************************************
% processOuput:  Process output and delete temporary files and variables.
% ------------------------------------------------------------------------------
% Called by: mads
% Calls:     getPollCenter, plotHistory
% VARIABLES:
%  Cache           = structure containing previously evaluated iterates
%    .Filter       =   structure containing filter data
%      .hmin       =     minimum allowable h value to be in Filter.F
%  BestF           = best feasible point found
%  BestI           = least infeasible point found
%  Problem         = structure containing optimization problem data
%    .File.P       =   name of parameter file
%    .nameCache    =   name of the base workspace Cache variable
%  Options         = structure of MADS parameters
%    .plotHistory1 =   flag for plotting f-value vs function evals
%    .hplothandle  =   handle of the history plot
%    .plotColor    =   string indicating color of history plot line
%  Param           = structure of output from user-provided Parameter file
%*******************************************************************************
function [Cache,BestF,BestI] = processOutput(Problem,Options)

% Close workspace environment
Cache = closeWorkspace(Problem);

% Plot performance history
if (Options.plotHistory1)
    plotHistory(Options.hplothandle,Options.plotColor,Cache);
end

% Store Best Solutions
BestF = getPollCenter(0,Cache.Filter,Cache);
BestI = getPollCenter(1,Cache.Filter,Cache);
if (BestF.h > Cache.Filter.hmin), BestF = []; end
if (BestI.h < Cache.Filter.hmin), BestI = []; end

% For any post-processing, call parameter file with BestF as an argument
if (exist(Problem.File.P,'file') == 2) && (nargin(Problem.File.P) < 0)
    Param = feval(Problem.File.P,BestF);
    setappdata(0,'PARAM',Param);
end
return

%*******************************************************************************
% closeWorkspace:  Shuts down base workspace appdata and deletes temp files.
% ------------------------------------------------------------------------------
% Called by: mads, ProcessOutput
% VARIABLES:
%  Cache           = structure containing previously evaluated iterates
%  Problem         = structure containing optimization problem data
%    .File.F       =   name of functions file
%    .nameCache    =   name of the base workspace Cache variable
%  surrName        = string used in constructing surrogate problem filename
%*******************************************************************************
function Cache = closeWorkspace(Problem)

% Store Cache and delete all base workspace appdata
Cache = getappdata(0,Problem.nameCache);
if isappdata(0,Problem.nameCache), rmappdata(0,Problem.nameCache); end
if isappdata(0,'SUR'),             rmappdata(0,'SUR');             end
if isappdata(0,'PARAM'),           rmappdata(0,'PARAM');           end
if isappdata(0,'P'),               rmappdata(0,'P');               end

% Delete function files used in search and surrogate optimization
lastwarn('');
if exist('gaPenalty','file'),  delete('gaPenalty.m'); end
if exist('SurObj','file'),     delete('SurObj.m'); end
if exist('SurNLConst','file'), delete('SurNLConst.m'); end
surrName = [Problem.File.F, '_Sur'];
if exist(surrName,'file'),           delete([surrName,'.m']); end
if exist([surrName,'_Cache'],'file'),delete([surrName,'_Cache.m']); end
if ~isempty(lastwarn)
    warning(['Search or surrogate file found but not deleted. \n', ...
            'File may have been moved from its original location. \n', ...
            'Please delete the file manually.'],[]);
end
return

%*******************************************************************************
% plotHistory:  Plot optimal value vs. number of function evaluations.
% ------------------------------------------------------------------------------
% Called by: update, ProcessOutput
% VARIABLES:
%  handle  = handle of the plot axes
%  color   = string code indicating color/linestyle of plot line
%  Cache   = Cache of previously computed iterates
%    .bfp  =   vector of best feasible f-values
%    .lip  =   vector of least infeasible f- and h-values
%    .size =   number of iterates in the Cache
%*******************************************************************************
function plotHistory(handle,color,Cache)

% Set up history plot axes handle
if isempty(handle)
    figure; handle = gca;
end

% Construct f plot using outside axes handle, if avaliable
axes(handle);
set([handle; get(handle,'Children')],'Visible','on');
if (any(isfinite(Cache.bfp)))
    plot(cumsum(Cache.nFunc(1:Cache.size),2), ...
        Cache.bfp(1:Cache.size)',  [color,'-'],'LineWidth',2);
else
    plot(cumsum(Cache.nFunc(1:Cache.size),2),...
        Cache.lip(1,1:Cache.size)',[color,'-'],'LineWidth',2);
end
title( 'Performance History', 'FontWeight','bold', 'FontSize',11);
xlabel('Number of Function Evaluations',           'FontSize',12);
ylabel('Objective Function Value  ',               'FontSize',12);
drawnow;
return;

%*******************************************************************************
% FUNCTIONS FOR STOCHASTIC PROBLEMS USING RANKING AND SELECTION (R&S)
%*******************************************************************************
%*******************************************************************************
% evalRSPointSet:  Evaluate a set of R&S candidate points
% ------------------------------------------------------------------------------
%     This function evaluates iterates to determine if feasible or infeasible.
%     If feasible, the feasible iterates are sent to be processed by RandSfeas;
%     and infeasible points are processes locally.
%     Code written by John Dunlap, Todd Sriver, and Mark Abramson
% ------------------------------------------------------------------------------
% Called by: evalPointSet
% Calls:     updateOmega, inOmega, getpID, isCacheHit, evalFunction, 
%            updateCache, updateFilter, RS
% VARIABLES:
%  P              = set of iterates with structure
%    .x           =   vector of continuous variables 
%    .p           =   cell array of categorical variables
%    .f           =   function value
%  rssuccess      = index of best trial point (0 = incumbent)
%  Filter         = updated filter
%    .feasible    =   indices of feasible iterates
%    .hmin        =   minimum h-value of an infeasible point
%  RSData         = structure of R&S parameters
%    .s0          =   number of initial replications per point
%    .iz          =   indifference zone parameter
%    .alpha       =   significance level alpha parameter
%    .nFuncLeft   =   number of function evaluations left before termination
%    .F           =   storage of function values in F
%  Problem        = structure containing problem parameters
%    .adaptOmega  =   flag for changing Omega of an MVp problem
%    .Omega       =   substructure that defines the feasible region
%    .isMVP       =   flag indicating the problem being solved is an MVP problem
%    .nameCache   =   name of the appdata holding the Cache
%    .File.F      =   name of the functions file
%    .fType       =   type of functions file (M=Matlab, F=Fotran, C=C/C++)
%  Pset           = set of iterates to be evaluated 
%  Options        = structure of user options
%  add2sur        = flag for including point in surrogate construction
%  feas           = vector of feasible iterates
%  infeas         = vector of infeasible iterates
%  nP             = number of points to be evaluated
%  pID            = unique ID for a specific categorical variable value
%  xnorm          = norm of the continuous variables
%  Cache          = structure of previously computed iterates
%    .isSurrPoint =   flags indicating if points are used to build surrogates
%    .size        =   number of iterates in the Cache
%  .nHits         =   number of Cache hits
%  hit            = flag indicating a Cache hit (previous computed iterate)
%  feasible       = flag indicating a point is feasible
%  fRSsuccess     = flag indicating R&S success with a feasible point
%  iRSsuccess     = flag indicating R&S success with an infeasible point
%  nf             = number of function evaluations
%  unfiltered     = flag indicating if current trial point is unfiltered
%  iterate        = temporary storage of current trial point
%*******************************************************************************
function [P,rssuccess,Filter,RSData] = ...
    evalRSPointSet(Problem,Pset,Options,Filter,RSData,add2sur)

% Initialize output variables
[P,feas,infeas] = deal([]);
rssuccess = 0;

% Evaluate each iterate, as appropriate
nP = length(Pset);
for k = 1:nP
    
    % update MVP Omega parameters, as is necessary
    if Problem.adaptOmega
        [Problem.Omega] = updateOmega(Problem,Options,Pset(k));
    end
    
    % Test if iterate lies in Omega
    if (inOmega(Pset(k),Problem.Omega))
        
        % Give trial point unique IDs to expedite Cache searching
        pID = 0;
        if Problem.isMVP
            pID = getpID(Pset(k).p,Problem.Omega.plist);
        end
        xnorm = norm(Pset(k).x,inf);
        
        % Test to see if the iterate is in the Cache
        Cache = getappdata(0,Problem.nameCache);
        hit = isCacheHit(Cache,Pset(k),xnorm,pID,50000);
        Cache.nHits = Cache.nHits + hit;
        setappdata(0,Problem.nameCache,Cache);
        
        % Evaluate function
        if (hit)
            P = [P,Pset(k)];
        else
            Pset(k) = evalFunction(Problem.File.F,Problem.fType,Pset(k), 0);
            feasible  = Pset(k).h <= Filter.hmin;    
            if (feasible)
                feas = [feas,Pset(k)];
            else
                infeas = [infeas,Pset(k)];
            end
        end
    end
end

[feas,fRSsuccess,Filter,RSData] = RS(Problem,feas,Filter,RSData);

% Update Cache and filter for feasible points
if ~isempty(feas)
    Cache = getappdata(0,Problem.nameCache);
    for j = 1:length(feas)
        pID = 0;
        if Problem.isMVP
            pID = getpID(feas(j).p,Problem.Omega.plist);
        end
        nf = feas(j).nFunc;
        Cache = updateCache(Cache,Problem, rmfield(feas(j),'nFunc'), ...
            norm(feas(j).x,inf),pID, nf);
        setappdata(0,Problem.nameCache,Cache);
        if (fRSsuccess == j)
            [unfiltered,Filter] = updateFilter(feas(j),Filter,Cache);
            if (add2sur == 2) || (add2sur == 1 && unfiltered)
                Cache = getappdata(0,Problem.nameCache);
                Cache.isSurrPoint(Cache.size) = 1;
                setappdata(0,Problem.nameCache,Cache);
            end
        end 
    end
    P = [P,rmfield(feas,'nFunc')];
end

% Evaluate and average infeasible points, update Cache and filter
iRSsuccess = 0;
unfiltered = 0;
if ~isempty(infeas)
    for i = 1:length(infeas)
        R(1,i) = infeas(i).f;
        for j = 2:RSData.s0
            iterate = evalFunction(Problem.File.F,Problem.fType,infeas(i), 0);
            R(j,i) = iterate.f;
        end
        pID = 0;
        if Problem.isMVP
            pID = getpID(infeas(i).p,Problem.Omega.plist);
        end
        nf = RSData.s0;
        infeas(i).f = mean(R(:,i));
        Cache = getappdata(0,Problem.nameCache);
        Cache = updateCache(Cache,Problem, infeas(i), ...
            norm(infeas(i).x,inf),pID, nf);
        setappdata(0,Problem.nameCache,Cache);
        if isfinite(infeas(i).h)
            [unfiltered,Filter] = updateFilter(infeas(i),Filter,Cache);
            if unfiltered, iRSsuccess = 1; end
            if (add2sur == 2) || (add2sur == 1 && unfiltered)
                Cache = getappdata(0,Problem.nameCache);
                Cache.isSurrPoint(Cache.size) = 1;
                setappdata(0,Problem.nameCache,Cache);
            end
        end
    end
    P = [P,infeas];
end
rssuccess = (fRSsuccess || iRSsuccess);

% Update R&S parameters
RSData.iz    = RSData.iz    * RSData.iz_rho;   
RSData.alpha = RSData.alpha * RSData.alpha_rho;
return;

%*******************************************************************************
% RS:  Execute the Screen and Selection Ranking & Selection (R&S) procedure.
% ------------------------------------------------------------------------------
%      See Pichitlamken, J. and B. L. Nelson, 
%      "Selection-of-the-best procedures for optimization via simulation",
%      Proceedings of the 2001 Winter Simulation Conference, pp. 401 - 407
%      Code written by John Dunlap, Todd Sriver,and Mark Abramson
% ------------------------------------------------------------------------------
% Called by: evalRSPointSet
% VARIABLES:
%  P              = set of iterates with structure
%    .x           =   vector of continuous variables 
%    .p           =   cell array of categorical variables
%    .f           =   function value
%  fRSsuccess     = index of best trial point (0 = incumbent)
%  Filter         = updated filter
%    .feasible    =   indices of feasible iterates
%  RSData         = structure of R&S parameters
%    .s0          =   number of initial replications per point
%    .iz          =   indifference zone parameter
%    .alpha       =   significance level alpha parameter
%    .nFuncLeft   =   number of function evaluations left before termination
%    .F           =   storage of function values in F
%  Problem        = structure containing optimization problem parameters
%  emptyFilter    = flag indicating if the current filter is empty
%  nP             = number of points to be evaluated
%  C              = temporary storage of the points to be evaluated
%  Cache          = Cache of data containing all previous computed iterates
%    .iterate     =   vector of iterates
%    .nFunc       =   numbers of function evaluations for each iterate
%  funcall        = number of function calls
%  F              = matrix of function values, by point and replication
%  nF             = number of points to be evaluated
%  S              = matrix of variances of the differences between the responses
%  lambda         = R&S parameter
%  a              = temporary storage
%  Nmax           = maximum number of replications required
%  meanvalue      = matrix of indices, mean responses, # function evaluations
%  ssmDone        = flag indicating termination criteria have been met
%  s              = replication counter
%  I              = indices of iterates still requiring more replications
%  Iold           = storage of previous iterations I matrix
%  R              = temporary storage of responses
%  minR           = minimum value in R
%*******************************************************************************
function [P,fRSsuccess,Filter,RSData] = RS(Problem,P,Filter,RSData)

fRSsuccess = 0;
if (isempty(P)), return; end
emptyFilter = isempty(Filter.feasible);

% Set up candidate points and store incumbent data
nP = length(P);
if emptyFilter
    [C(1:nP).x] = deal(P(1:nP).x);
    [C(1:nP).p] = deal(P(1:nP).p);
    [C(1:nP).f] = deal(P(1:nP).f);
else
    Cache = getappdata(0,Problem.nameCache);
    C(1).x = Cache.iterate(Filter.feasible(1)).x;
    C(1).p = Cache.iterate(Filter.feasible(1)).p;
    C(1).f = Cache.iterate(Filter.feasible(1)).f;
    [C(2:nP+1).x] = deal(P(1:nP).x);
    [C(2:nP+1).p] = deal(P(1:nP).p);
    [C(2:nP+1).f] = deal(P(1:nP).f);
end
funcall = nP;

%------------- STEP 1: Initialization ---------------

% Initially evaluate each candidate point
F(1,:) = [C.f];
nF = size(F,2);
for i = 1:nF
    if (i == 1) && (size(RSData.F,1) >= RSData.s0)
        F(1:RSData.s0,1) = RSData.F(1:RSData.s0,1);
    else
        for j = 2:RSData.s0
            iterate = evalFunction(Problem.File.F,Problem.fType,C(i), 0);
            F(j,i)  = iterate.f;
            funcall = funcall + 1;
        end
    end
end

% Store in S the variance for each candidate pair
S = zeros(nF);
for i = 1:nF
    for j = 1:i-1
        S(i,j) = var(F(:,i)-F(:,j));
        S(j,i) = S(i,j);
    end
end

% ------------ STEP 2: Procedure parameters --------------

% Compute maximum number of replications required
lambda = RSData.iz / 2;            %Lambda parameter recommended on pp. 403
a = ((RSData.s0-1)*S / (4*(RSData.iz-lambda)))* ...
    (((nF-1)/(2*(RSData.alpha)))^(2/(RSData.s0-1))-1);
a = a - diag(diag(a));
Nmax = max(max(floor(a/lambda)));

% Set termination flag if initial sample was enough; otherwise, next stage
meanvalue = [];           %Matrix to store indices, mean responses, & nFunc         
if (RSData.s0 > Nmax)
    meanvalue = [meanvalue; (1:nF)', mean(F,1)', RSData.s0*ones(nF,1)];
    ssmDone = 1;
else
    s = RSData.s0;
    I = (1:nF)';             % stores set of candidates still in play
    ssmDone = 0;
end

% If reached max function calls
if ((funcall >= RSData.nFuncLeft) && ~ssmDone)
    meanvalue = [meanvalue; I, mean(F(:,I))', s*ones(nF,1)];
    ssmDone = 1;
end

% Begin main loop
while (~ssmDone)
    
    % ------------ STEP 3: Screening --------------
    
    Iold = I;
    I = []; 
    R = zeros(nF,1);
    R(Iold) = R(Iold) + sum(F(:,Iold))';
    
    for i = 1:length(Iold)
        minR = realmax;
        if (R(Iold(i))==Inf)
            meanvalue = [meanvalue; Iold(i), mean(F(:,Iold(i))), s];
        else
            for j = 1:length(Iold)
                if (i ~= j)
                    minR = min(minR, R(Iold(j)) + a(Iold(i),Iold(j)));
                end
            end
            if (R(Iold(i)) <= minR - s*lambda) % Note:<= because minimization
                I = [I ; Iold(i)];
            else
                meanvalue = [meanvalue; Iold(i), mean(F(:,Iold(i))), s];
            end
        end
    end
    
    % ------------ STEP 4: Stopping Rule --------------
    
    if isempty(I)                           % I is empty
        ssmDone = 1;
    elseif (length(I) == 1)                 % Only one candidate left
        ssmDone = 1;
        meanvalue = [meanvalue; I(1), mean(F(:,I(1))), s];
        
        % Take new observation from screening survivors, check extra stopping rule
    else
        s = s + 1;
        for i = 1:length(I)
            if (I(i) == 1) && (size(RSData.F,1) >= s)
                F(s,1) = RSData.F(s,1);
            else
                iterate = evalFunction(Problem.File.F,Problem.fType,C(i), 0);
                F(s,I(i)) = iterate.f;
                funcall = funcall + 1;
            end
        end
        
        % Terminate: enough samples have been taken or max function calls exceeded
        if (s == Nmax + 1) || (funcall >= RSData.nFuncLeft)
            meanvalue = [meanvalue; I, mean(F(:,I))', s*ones(length(I),1)];
            ssmDone = 1;
        end
    end
end

% Sort meanvalue matrix by mean
meanvalue = sortrows(meanvalue,2);

if emptyFilter
    for i = 1:size(meanvalue,1)
        P(meanvalue(i,1)).f     = meanvalue(i,2);
        P(meanvalue(i,1)).nFunc = meanvalue(i,3);
    end
    fRSsuccess = meanvalue(1,1);          % Stores the index of best candidate
else
    if (meanvalue(1,1) ~= 1)
        fRSsuccess = (meanvalue(1,1)-1);
    end
    
    % Store the average function value and number of function evaluations
    for i = 1:size(meanvalue,1)
        if (meanvalue(i,1) ~= 1)
            P(meanvalue(i,1)-1).f     = meanvalue(i,2);
            P(meanvalue(i,1)-1).nFunc = meanvalue(i,3);
        end
    end
end

% Replace incumbent
if emptyFilter
    RSData.F = F(1:meanvalue(1,3),meanvalue(1,1));  
    
    %Replaces current history based on the number of replications for new
    %incumbent.  Since it is not changed outside this function, it must be
    %updated here.
elseif (meanvalue(1,1) ~= 1) 
    RSData.F = F(1:meanvalue(1,3),meanvalue(1,1));  
    if (meanvalue(1,1) ~= 1)
        for i = 1:size(meanvalue,1)
            if (meanvalue(i,1) == 1)
                nf = meanvalue(i,3);
                temp = meanvalue(i,2);
            end
        end
        if (nf > Cache.nFunc(Filter.feasible(1)))      %Only replace if more
            Cache = getappdata(0,Problem.nameCache);    %observations
            Cache.nFunc(Filter.feasible(1)) = nf;
            Cache.iterate(Filter.feasible(1)).f = temp;
            setappdata(0,Problem.nameCache,Cache);
        end
    end
    
    % If more data is available for incumbent
elseif (size(RSData.F,1) < size(F(:,meanvalue(1,1)),1))
    Cache = getappdata(0,Problem.nameCache);
    Cache.iterate(Filter.feasible(1)).f = meanvalue(1,2);
    Cache.nFunc(Filter.feasible(1)) = meanvalue(1,3);
    setappdata(0,Problem.nameCache,Cache);
    RSData.F = F(1:meanvalue(1,3),meanvalue(1,1));
end
return

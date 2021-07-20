%==========================================================================
% THE USER CAN DEFINE THE PROBLEM AND SET OPTIONS IN THE FOLLOWING LINES:                                                       
%==========================================================================

function [modelname,paths,opts,submodels,prev_ident_pars] = options()

global ar

%%% (1) MODEL: 
modelname = ar.ia.modelname;

%%% (2) PATHS:
paths.meigo     = ar.ia.paths.meigo;      
paths.models    = ar.ia.paths.models;
paths.results   = ar.ia.paths.results;
paths.functions = ar.ia.paths.functions;
                            
%%% (3) IDENTIFIABILITY OPTIONS:
opts.numeric    = ar.ia.opts.numeric;          % calculate rank numerically (= 1) or symbolically (= 0)
opts.replaceICs = ar.ia.opts.replaceICs;       % replace states with specific initial conditions (= 1) or use generic values (= 0) when calculating rank
opts.checkObser = ar.ia.opts.checkObser;       % check state observability, i.e. identifiability of initial conditions (1 = yes; 0 = no).
opts.checkObsIn = ar.ia.opts.checkObsIn;       % check input observability (1 = yes; 0 = no).
opts.unidentif  = ar.ia.opts.unidentif;        % use method to try to establish unidentifiability instead of identifiability, when using decomposition. 
opts.forcedecomp= ar.ia.opts.forcedecomp;      % always decompose model (1 = yes; 0 = no).
opts.decomp     = ar.ia.opts.decomp;           % decompose model if the whole model is too large (1 = yes; 0 = no: instead, calculate rank with few Lie derivatives).
opts.decomp_user= ar.ia.opts.decomp_user;      % when decomposing model, use submodels specified by the user (= 1) or found by optimization (= 0). 
opts.maxLietime = ar.ia.opts.maxLietime;       % max. time allowed for calculating 1 Lie derivative.
opts.maxOpttime = ar.ia.opts.maxOpttime;       % max. time allowed for every optimization (if optimization-based decomposition is used).
opts.maxstates  = ar.ia.opts.maxstates;        % max. number of states in the submodels (if optimization-based decomposition is used).
opts.nnzDerU    = ar.ia.opts.nnzDerU;          % numbers of nonzero derivatives of the measured inputs (u); may be 'inf'
opts.nnzDerW    = ar.ia.opts.nnzDerW;          % numbers of nonzero derivatives of the unmeasured inputs (w); may be 'inf'
opts.nnzDerIn   = ar.ia.opts.nnzDerIn;         % deprecated option

%%% (4) AFFINE OPTIONS (to use the ORC-DF algorithm):
opts.affine              = ar.ia.opts.affine;             % use the ORC-DF algorithm for affine control systems (=1) or not(=0).
opts.affine_tStage       = ar.ia.opts.affine_tStage;      % max. computation time for the last iteration.
opts.affine_kmax         = ar.ia.opts.affine_kmax;        % max. number of iterations.
opts.affine_parallel     = ar.ia.opts.affine_parallel;    % use parallel toolbox (=1) or not (=0) to calculate partial ranks.
opts.affine_workers      = ar.ia.opts.affine_workers;     % number of workers for parallel pool.
opts.affine_graphics     = ar.ia.opts.affine_graphics;    % display graphics (=1) or not (=0)
opts.affine_delete_model = ar.ia.opts.affine_delete_model;% delete affine model when finished (=1) or not (=0).

%%% (5) DECOMPOSITION OPTIONS -- User-defined submodels for decomposition (may be left = []): 
submodels = ar.ia.opts.submodels; 
%- Submodels are specified as a vector of states, as e.g.:
%         submodels{1}  = [1 2];
%         submodels{2}  = [1 3];

%%% (6) MULTI-EXPERIMENT OPTIONS:
opts.multiexp              = ar.ia.opts.multiexp;               % Execute multi-experiment analysis (=1) or not (=0).
opts.multiexp_numexp       = ar.ia.opts.multiexp_numexp;        % Number of experiments.
opts.multiexp_user_nnzDerU = ar.ia.opts.multiexp_user_nnzDerU;  % Set manually the number of non-zero known input derivatives in each experiment (=1) or not (=0).
opts.multiexp_nnzDerU      = ar.ia.opts.multiexp_nnzDerU;       % Number of non-zero known input derivatives in each experiment (Rows=inputs;Columns=experiments).
opts.multiexp_user_nnzDerW = ar.ia.opts.multiexp_user_nnzDerW;  % Set manually the number of non-zero unknown input derivatives in each experiment (=1) or not (=0).
opts.multiexp_nnzDerW      = ar.ia.opts.multiexp_nnzDerW;       % Number of non-zero unknown input derivatives in each experiment (Rows=inputs;Columns=experiments).
opts.multiexp_user_ics     = ar.ia.opts.multiexp_user_ics;      % Set manually the initial conditions for each experiment (=1) or not (=0).
%- Multi-experiment initial conditions (Rows=variables;Columns:experiments):
opts.multiexp_ics       = ar.ia.opts.multiexp_ics;
%- Which initial conditions are replaced (Rows=variables;Columns=experiments):
opts.multiexp_known_ics = ar.ia.opts.multiexp_known_ics;

%%% (7) LIE SYMMETRIES & REPARAMETERIZATION OPTIONS:
opts.ansatz     = ar.ia.opts.ansatz;     % Type of Ansatz: %  uni -> Univariate (1)
                                         %  par -> Partially variate (2)
                                         %  multi -> Multivariate (3)
opts.degree     = ar.ia.opts.degree;     % Degree of Ansatz Polynomial
opts.tmax       = ar.ia.opts.tmax;       % Maximum degree of Lie series
opts.ode_n      = ar.ia.opts.ode_n;      % (Only in Matlab R2020a and later:) use ode solver (=1) or not (=0)
opts.use_existing_results = ar.ia.opts.use_existing_results;  % if the model has already been analysed with STRIKE-GOLDD (=1) or not (=0)
opts.results_file = ar.ia.opts.results_file;                  % .mat file to use if use_existing_results = 1                              

%%% (8) KNOWN/IDENTIFIABLE PARAMETERS (parameters assumed known, or already classified as identifiable):
prev_ident_pars = ar.ia.opts.prev_ident_pars;
%- example:
% syms p2 p5
% prev_ident_pars = [p2 p5];


end

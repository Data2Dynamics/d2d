%==========================================================================
% Matlab toolbox for structural identifiability and observability
% analysis of nonlinear models
%--------------------------------------------------------------------------
% based on:
% StrIkE-GOLDD (v3.0, last modified: 19/10/2020)
% https://github.com/afvillaverde/strike-goldd
%--------------------------------------------------------------------------
% options can be specified in 'ar.ia.opts'



% m     model number in ar struct. default value is m=1

function arSIAInit(m)

global ar

if isfield(ar,'ia')
    if isfield(ar.ia,'opts')
        opts_old = ar.ia.opts;
        opts_old.use_existing_results = 0;
    end
    ar = rmfield(ar,'ia');
end

%userManual = [ar.ia.paths.strike_goldd '/doc/STRIKE-GOLDD_manual.pdf'];

disp(' ')
disp('====================================================================')
disp('    this toolbox uses <a href="https://github.com/afvillaverde/strike-goldd">STRIKE-GOLDD (v3.0)</a> package for analysing')
disp('the structural identifiability and observability of nonlinear models')
%disp('                         (<a href="matlab: open(userManual)">user manual</a>)')
disp('--------------------------------------------------------------------')
disp('   For more information about the theory of finding and breaking');
disp('                   Lie Symmetries click <a href="https://doi.org/10.3390/sym12030469">here</a>')
disp('====================================================================')


if ~exist('ar','var') && ~isfield(ar,'model')
    error('Please initialize and load the model');
end

if ~exist('m','var') || isempty(m)
    m = 1;
end

if isempty(ar.model(m).y)
    error('The toolbox is not compatible with models with no observables!');
end

if exist('ar','var') && isfield(ar,'ia') && isfield(ar.ia,'paths')
    cd(ar.ia.paths.d2dmodel)
end
    

% create the model
arStrikeGolddModel(m);

% set options
if exist('opts_old','var')
    ar.ia.opts = opts_old;
else
    InitiateOptions;
end


ar.ia.about = 'this field contains information used in the strike-goldd toolbox';


currentFolder = pwd;
ar.ia.paths.d2dmodel = currentFolder;


ar_path = fileparts(which('arInit.m'));
strike_goldd_path = strcat(ar_path,filesep,'ThirdParty',filesep,'strike-goldd',filesep,'STRIKE-GOLDD');
ar.ia.paths.strike_goldd = strike_goldd_path;

ID = 1;      % check if the model is already exist
if ~exist('Results', 'dir')
    ID = 0;
    status = mkdir('Results');
    if status == 0
        error('Making folder "Results" was not successful');
    end
end

model_path = strcat(currentFolder,filesep,'Results');
cd(model_path)

if ~exist('Identifiability_Analysis', 'dir')
    ID = 0;
    status = mkdir('Identifiability_Analysis');
    if status == 0
        error('Making folder "Identifiability_Analysis" was not successful');
    end
end
model_path = strcat(model_path,filesep,'Identifiability_Analysis');
cd(model_path)

if ~exist(ar.ia.modelname_checkstr, 'dir')
    ID = 0;
    status = mkdir(ar.ia.modelname_checkstr);
    if status == 0
        error('Making folder "%s" was not successful',ar.ia.modelname_checkstr);
    end
end
model_path = strcat(model_path,filesep,ar.ia.modelname_checkstr);
cd(model_path)

if ~exist('result', 'dir')
    ID = 0;
    status = mkdir('result');
    if status == 0
        error('Making folder "%s" was not successful','result');
    end
end

if ~exist(strcat(ar.ia.modelname,'.m'),'file') || ~exist(strcat(ar.ia.modelname,'.mat'),'file')
    ID = 0;   
    % write the model to file
    arStrikeGolddWriteModel   
    % run the model
    run([ar.ia.modelname,'.m'])
end

% check if the result is already exist
if ID==1
    folderInfo={};
    tmp = dir('result');
    for i=1:length(tmp)
        folderInfo{end+1} = tmp(i).name;
    end
    status = 0;
    for i=1:length(folderInfo)
        status = contains(folderInfo,strcat('id_results_',ar.ia.modelname));
    end
    if sum(status)==1
        ar.ia.opts.use_existing_results = 1;
        ar.ia.opts.results_file = folderInfo{status};
        fprintf(' * The model has already been analysed!\n\n');
    end
end


% redirect to the strike_goldd path
cd(strike_goldd_path)

% create needed paths
ar.ia.paths.meigo     = strcat(pwd,filesep,'MEIGO64/MEIGO');
ar.ia.paths.models    = strcat(model_path);
ar.ia.paths.results   = strcat(model_path,filesep,'result');
ar.ia.paths.functions = strcat(pwd,filesep,'functions');

% initiate strike-goldd and adds paths to the matlab search path
install

% add paths to the matlab search path
addpath(model_path)
addpath(fullfile(model_path,'result'))
ar.ia.addedpaths{end+1} = model_path;
ar.ia.addedpaths{end+1} = fullfile(model_path,'result');

ar.ia = orderfields(ar.ia);

% set the options
options;

% remove paths from the matlab search path.
rmpath(ar.ia.addedpaths{:})

% redirect to the original folder
cd(ar.ia.paths.d2dmodel)

fprintf(' Initialization ... [done]\n');

end



function InitiateOptions
% defult values to uses STRIKE-GOLDD (v3.0)
% these values can be changes by the user after initialization.
% this function initializes all the information needed in option.m

global ar

% AUTOMATIC REPARAMETERIZATION
ar.ia.opts.autoRepar = 0;

% IDENTIFIABILITY OPTIONS:
ar.ia.opts.numeric     = 0;       % calculate rank numerically ( = 1) or symbolically ( = 0)
ar.ia.opts.replaceICs  = 0;       % replace states with specific initial conditions ( = 1) or use generic values ( = 0) when calculating rank
ar.ia.opts.checkObser  = 1;       % check state observability, i.e. identifiability of initial conditions (1  = yes; 0  = no).
ar.ia.opts.checkObsIn  = 1;       % check input observability (1  = yes; 0  = no).
ar.ia.opts.unidentif   = 0;       % use method to try to establish unidentifiability instead of identifiability, when using decomposition.
ar.ia.opts.forcedecomp = 0;       % always decompose model (1  = yes; 0  = no).
ar.ia.opts.decomp      = 0;       % decompose model if the whole model is too large (1  = yes; 0  = no: instead, calculate rank with few Lie derivatives).
ar.ia.opts.decomp_user = 0;       % when decomposing model, use submodels specified by the user ( = 1) or found by optimization ( = 0).
ar.ia.opts.maxLietime  = 100;     % max. time allowed for calculating 1 Lie derivative.
ar.ia.opts.maxOpttime  = 30;      % max. time allowed for every optimization (if optimization-based decomposition is used).
ar.ia.opts.maxstates   = 6;       % max. number of states in the submodels (if optimization-based decomposition is used).
ar.ia.opts.nnzDerU     = inf;     % numbers of nonzero derivatives of the measured inputs (u); may be 'inf'
ar.ia.opts.nnzDerW     = 1;       % numbers of nonzero derivatives of the unmeasured inputs (w); may be 'inf'
ar.ia.opts.nnzDerIn    = ar.ia.opts.nnzDerU; % deprecated option

% AFFINE OPTIONS (to use the ORC-DF algorithm):
ar.ia.opts.affine               = 0;     % use the ORC-DF algorithm for affine control systems ( =1) or not( =0).
ar.ia.opts.affine_tStage        = 1000;  % max. computation time for the last iteration.
ar.ia.opts.affine_kmax          = 4;     % max. number of iterations.
ar.ia.opts.affine_parallel      = 0;     % use parallel toolbox ( =1) or not ( =0) to calculate partial ranks.
ar.ia.opts.affine_workers       = 4;     % number of workers for parallel pool.
ar.ia.opts.affine_graphics      = 1;     % display graphics ( =1) or not ( =0)
ar.ia.opts.affine_delete_model  = 1;     % delete affine model when finished ( =1) or not ( =0).

% DECOMPOSITION OPTIONS -- User-defined submodels for decomposition (may be left  = []):
ar.ia.opts.submodels  = [];
%- Submodels are specified as a vector of states, as e.g.:
%         submodels{1}   = [1 2];
%         submodels{2}   = [1 3];

% MULTI-EXPERIMENT OPTIONS:
ar.ia.opts.multiexp               = 0;     % Execute multi-experiment analysis ( =1) or not ( =0).
ar.ia.opts.multiexp_numexp        = 2;     % Number of experiments.
ar.ia.opts.multiexp_user_nnzDerU  = 0;     % Set manually the number of non-zero known input derivatives in each experiment ( =1) or not ( =0).
ar.ia.opts.multiexp_nnzDerU       = [0 1]; % Number of non-zero known input derivatives in each experiment (Rows =inputs;Columns =experiments).
ar.ia.opts.multiexp_user_nnzDerW  = 0;     % Set manually the number of non-zero unknown input derivatives in each experiment ( =1) or not ( =0).
ar.ia.opts.multiexp_nnzDerW       = [1 0]; % Number of non-zero unknown input derivatives in each experiment (Rows =inputs;Columns =experiments).
ar.ia.opts.multiexp_user_ics      = 0;     % Set manually the initial conditions for each experiment ( =1) or not ( =0).
%- Multi-experiment initial conditions (Rows =variables;Columns:experiments):
ar.ia.opts.multiexp_ics        = [ [1,0,0,1,0,1,0,0,0].', [1,0,0,1,0,1,0,0,0].' ];
%- Which initial conditions are replaced (Rows =variables;Columns =experiments):
ar.ia.opts.multiexp_known_ics  = [ [0,1,1,0,1,0,1,1,1].', [0,1,1,0,1,0,1,1,1].' ];

% LIE SYMMETRIES & REPARAMETERIZATION OPTIONS:
ar.ia.opts.ansatz      = 2; % Type of Ansatz: %  uni -> Univariate (1)
%  par -> Partially variate (2)
%  multi -> Multivariate (3)
ar.ia.opts.degree      = 2; % Degree of Ansatz Polynomial
ar.ia.opts.tmax        = 4; % Maximum degree of Lie series
ar.ia.opts.ode_n       = 1; % (Only in Matlab R2020a and later:) use ode solver ( =1) or not ( =0)
ar.ia.opts.use_existing_results  = 0; % if the model has already been analysed with STRIKE-GOLDD ( =1) or not ( =0)
%ar.ia.opts.results_file  = 'id_results_1D_BIG_p_16-Oct-2020.mat'; % .mat file to use if use_existing_results  = 1
ar.ia.opts.results_file  = ''; % .mat file to use if use_existing_results  = 1

% KNOWN/IDENTIFIABLE PARAMETERS (parameters assumed known, or already classified as identifiable):
ar.ia.opts.prev_ident_pars  = [];
%- example:
% syms p2 p5
% prev_ident_pars = [p2 p5];

end


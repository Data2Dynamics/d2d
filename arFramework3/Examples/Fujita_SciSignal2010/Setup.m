%% Compile model
arInit

warning off; % switch off warnings -> caution in equilibrium equations
arLoadModel('model_AktPathwayFujita');
warning on;  % switch on warnings

% Load experimental data
for n = 1:6
    arLoadData(['experimentaldata' num2str(n) ],1,'csv');
end

% Compile model
arCompileAll;
arFindInputs;
%% Constraint parameters
ar.lb = -8 * ones(size(ar.lb)); % lower parameter bounds
ar.ub =  8 * ones(size(ar.ub)); % upper parameter bounds

%% Numerical settings
ar.config.atol     = 1e-5;
ar.config.rtol     = 1e-5;
ar.config.maxsteps = 1e5;

% Optimizer settings
ar.config.optim.TolFun           = 1e-5;
ar.config.optim.PrecondBandWidth = inf;
ar.config.optim.Display          = 'iter';
ar.config.optim.MaxIter          = 1e4;

arLoadPars('ParamsFujita2010')
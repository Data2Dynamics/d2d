% Shaping the dose response of IL6 induced target gene transcription
addpath('Helper');
arInit;
arLoadModel('il6_model');

% Load one actual dataset to make sure we loaded the right parameters
arLoadData('bohl/hep_2005_02_03_Cont40ng_T90min', 1, [], true);

%% predictions
arLoadData('braun/prediction/APP_variable_prediction_TC');
arLoadData('steadystate/steadystate',1,'csv',true)

%% compile
arCompileAll(true);

% Tight equilibration tolerances
ar.config.eq_tol=1e-10;
ar.config.atol = 1e-8;
ar.config.rtol = 1e-8;
dataSets;

% Pre-equilibration settings (use a steady state condition to equilibrate
% the system)
arClearEvents;
arFindInputs;
arSteadyState(1,arFindCondition(ar,'steady'),'all');

%% settings
ar.config.optimizer = 5;
ar.config.optim.TolX = 1e-8;
ar.config.optim.MaxIter = 1000;
ar.config.optim.Display = 'iter';
ar.config.maxsteps = 1000;
ar.config.fiterrors = 1;
arLoadPars( 'finalized_model' );

fn = 'il6_experiment_design';
save(fn, 'ar');

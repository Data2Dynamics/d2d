%% initialize full model
arInit
arLoadModel('FullModel');
arLoadData('TimeCourseData', 1);
arCompileAll;

% set parameters for parameter estimation and optimization
ar.lb(:)   = -5;
ar.ub(:)   = 5;
ar.config.atol = 1e-8;
ar.config.rtol = 1e-8;
ar.config.optim.TolFun = 1e-8;
ar.config.optim.TolX = 1e-8;
ar.config.maxsteps = 1e5;

% set antibody specificity and STAT5A/B ratio
arSetPars('specC17',0.107, 0, 0);
arSetPars('ratio',0.693, 0, 0); 

% load parameters from best fit of LHS(500)
arLoadPars('FullModel')


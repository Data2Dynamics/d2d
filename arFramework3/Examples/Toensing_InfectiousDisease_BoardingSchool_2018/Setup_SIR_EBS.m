%initialize framework, load model and data 
arInit;
arLoadModel('SIR');
arLoadData('English_Boarding_School_1978');

% compile model and data 
arCompileAll;

% Extend upper parameter boundaries (i.e. parameter search space)
ar.ub(:)=4;

% set optimizer tolerances
arSetOptimTol;
ar.config.optim.MaxFunEvals = 5000;
ar.config.optim.MaxIter = 5000;

% set integrator tolerances
ar.config.rtol = 1e-7;
ar.config.atol = 1e-7;
 
arLoadPars('BestFit_EBS');


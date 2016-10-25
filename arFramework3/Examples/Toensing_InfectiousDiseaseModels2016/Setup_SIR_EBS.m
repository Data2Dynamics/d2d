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
 
% sinlge fit
arFit;

% multistart fit
arFitLHS(100,123);
% or load best fit parameters
% arLoadPars('BestFit_EBS');

% show parameter values
arPrint;

% plot results of multistart fit
arPlotFits

% plot observables, internal sates and fluxes
arQplot('xyv',1);
arPlot;

% save results
arSave('SIR_EBS_Multistart100')

% likelihood profiles
arPLEInit

% Set tolerances
pleGlobals.relchi2stepincrease(5) = 0.01;
pleGlobals.minstepsize(:) = 1e-4;

% calculate profiles
ple(1:5,200)

% plot profiles
plePlotMulti;

% plot trajectories along profiles
arPLETrajectories;

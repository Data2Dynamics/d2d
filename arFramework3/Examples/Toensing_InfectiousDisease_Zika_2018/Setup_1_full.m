%initialize framework, load model and data 
arInit;
arLoadModel('Zika_full');
arLoadData('Zika_Colombia');

% compile model and data 
arCompileAll;

% set parameter ranges
arSetParsPattern('init_',[],1,1,-5,10)
arSetParsPattern('beta_',[],1,1,-10,5)
arSetParsPattern('nu_',[],1,1,-10,5)
arSetParsPattern('gamma_',[],1,1,-10,5)
arSetPars('mu_v',[],1,1,-10,5)
arSetPars('sd_rel',-1,1,1,-5,-0.3)
arSetPars('sd_abs',-1,1,1,-5,2)
arSetPars('kappa_as',0.8,1,0,0,1)

% set optimizer tolerances
arSetOptimTol;
ar.config.optim.MaxIter = 5000;
ar.config.optim.MaxFunEvals = 5000;
ar.config.atol  =  1.0000e-09;
ar.config.rtol  =  1.0000e-09;
ar.config.maxsteps = 5000;
ar.config.atolV = 1;
ar.config.atolV_Sens = 1;


% sinlge fit
arFit;

%% The following code indicates Multi-Start fitting, Profile likelihood calculation and Uncertainty analysis
% % multistart fit
% arFitLHS(1000,1122);
% % or load best fit parameters
arLoadPars('BestFit_zika_full');
% 
% % show parameter values
% arPrint;
% 
% % plot results of multistart fit
% arPlotFits
% 
% % plot observables, internal sates and fluxes
% arQplot('xyv',1);
% arPlot;
% 
% % save results
% arSave('SIR_Zika_full_Multistart1000')
% 
% % likelihood profiles
% arPLEInit
% 
% % Set tolerances
% ar.ple.minstepsize(6) = 1e-5;
% 
% % calculate profiles
% ple(find(ar.qFit==1),1e3);
% 
% % extend profile suntin they hit the bound
% pleExtend;
% 
% % plot profiles
% plePlotMulti;
% 
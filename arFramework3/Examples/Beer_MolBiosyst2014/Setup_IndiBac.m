close all; clc;

% initialize model
arInit
arLoadModel('TwoStepSynthesisWithGlu');
arLoadData('20131208_indC_native_sfp_rep4_siData_1', 1,'csv');
arLoadData('20131208_synT1_sfp_siData_1', 1,'csv');
arLoadData('20131208_synT3_sfp_siData_1', 1,'csv');
arLoadData('20131208_synT5_sfp_siData_1', 1,'csv');
arCompileAll;

% set parameters for parameter estimation and optimization
ar.lb(:)   = -5;
ar.ub(:)   = 5;
ar.config.atol = 1e-8;
ar.config.rtol = 1e-8;
ar.config.optim.TolFun = 1e-8;
ar.config.optim.TolX = 1e-8;
ar.config.maxsteps = 1e5;

% load parameters from best fit
arLoadPars('BestFit')

% plot data sets shown in publication
ar.model.qPlotYs([1, 2, 8, 14]) = 1; 

arPlot

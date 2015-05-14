%% Full model of Boehm et al. Journal of proteome research 2014

close all; clear all; clc;

% initialize full model
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

arPlot



% phase plane plot

figure

% experimental data
pADat = ar.model.data.yExp(:,1);
pBDat = ar.model.data.yExp(:,2);
rADat = ar.model.data.yExp(:,3);
ptotDat = (rADat .* pADat + ((100 - rADat) .* pBDat)) / 100;

% full model simulations
pAFull = ar.model.data.yFineSimu(:,1);
pBFull = ar.model.data.yFineSimu(:,2);
rAFull = ar.model.data.yFineSimu(:,3);
ptotFull = (rAFull .* pAFull + ((100 - rAFull) .* pBFull)) / 100;

% plot
plot(ptotDat, rADat, 'ok');
xlabel('STAT5A phosphorylation degree [%]');
ylabel('relative STAT5A amount [%]');
axis([0 100 0 100])
hold on
plot(ptotFull, rAFull, 'k');
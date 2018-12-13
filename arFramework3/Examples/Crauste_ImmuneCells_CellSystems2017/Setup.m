close all
clc

% Load models & data
arInit;
arLoadModel('ModelImmuneCells');
arLoadData('DataCrausteImmune');
arCompileAll;

ar.config.fiterrors = -1;

%Set fmincon interior-point optimizer
ar.config.maxsteps = 4.e4;
ar.config.optimizer = 1;
ar.config.atol = 1.e-10;
ar.config.rtol = 1.e-10;
ar.config.optim.MaxIter = 1.e4;
ar.config.optim.InitTrustRegionRadius = 1;

arLoadPars('bestFit')
arPrint
arPlot
% Load models & data
clear all;
arInit

arLoadModel('ABC_model');
arLoadData('ABC_data_BCobs'); %Data with equidistant observation
%of state B and C for t=0,10,..100
%arLoadData('ABC_data_B_sparseObs'); %Data with sparse observation of
                                 %state B
arCompileAll();

%!Set laplace distribution!
ar.config.optimizer = 19;
ar.config.optimizers = [ar.config.optimizers 'fmincon_laplace'];

%Fit without estimated errors
%Take error 0.1 of data def
ar.config.fiterrors = -1;
arSetPars('sd_B_au',[],0);
arSetPars('sd_C_au',[],0);
arFitLHS(100);
arPlotFits
 
arSave('ABC_laplace')
 
arPLEInit
ple(1:3)

%% Second part

%Fit with estimated errors
ar.config.fiterrors=0;
ar.model(1).data(1).yExpStd(:,2) = NaN;
arSetPars('sd_B_au',[],1);
arFitLHS(100)
arPlotFits
arSave('ABC_laplace_errorEstimate')
arPLEInit
ple(1:4)
% 
% %If you want to check the error distribution
% %check_laplace
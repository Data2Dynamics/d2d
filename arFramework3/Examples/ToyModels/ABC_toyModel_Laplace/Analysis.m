% Load models & data
clear all;
Setup

%%
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
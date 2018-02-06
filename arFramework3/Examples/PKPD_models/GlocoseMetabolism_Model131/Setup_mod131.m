%% source: DDmore (http://repository.ddmore.foundation/model/DDMODEL00000131#Overview)

close all; clc;

% initialize model
arInit
arLoadModel('mod131');
arLoadData('SimuInsulinAction','mod131')
arCompileAll;


% initial values from DDmore (ERROR_prop is zero)
ar.qLog10(3) = 0;
ar.qFit(3) = 0;
varInitWOinit = [log10(150),log10(0.004),0,log10(0.03),log10(0.3),log10(0.03),log10(1),log10(0.3),log10(0.02),log10(0.4),log10(0.5)];
ar.p = varInitWOinit;


arFitLHS(1000)
arPlotFits

arPlot
hold on
plot([245,245],[0.2,0],':k')    % mark the start of the first insulin infusions
plot([145,145],[0.2,0],':k')    % mark the start of the second insulin infusions
arPrint


%% Profile-Likelihood
arPLEInit
ple
plePlotMulti

%% source: DDmore (http://repository.ddmore.foundation/model/DDMODEL00000131#Overview)

close all; clc;

Setup_mod131
%%
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

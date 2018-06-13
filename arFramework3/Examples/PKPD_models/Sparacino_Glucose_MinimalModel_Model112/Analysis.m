%% source: DDmore (http://repository.ddmore.foundation/model/DDMODEL00000112#Overview)

close all; clc;

Setup_mod112

%% fit
arFitLHS(50)

arPlotFits
arPlot
arPrint



%% profile likelhood
arPLEInit
ple
plePlotMulti


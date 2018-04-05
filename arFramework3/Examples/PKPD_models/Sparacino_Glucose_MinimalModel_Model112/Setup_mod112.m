%% source: DDmore (http://repository.ddmore.foundation/model/DDMODEL00000112#Overview)

close all; clc;

% initialize model
arInit
arLoadModel('mod112');
arLoadData('data112');
arCompileAll;


% fit
arFitLHS(50)

arPlotFits
arPlot
arPrint



%% profile likelhood
arPLEInit
ple
plePlotMulti


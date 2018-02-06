%% source: DDmore (http://repository.ddmore.foundation/model/DDMODEL00000132#Overview)

close all; clc;

% initialize model
arInit
arLoadModel('mod132');
arLoadData('Real_minimalModel')
arCompileAll;


arFitLHS(100)
arPlotFits

arPlot
arPrint

%% Profile-Likelhood for all parameters
arPLEInit
ple
plePlotMulti



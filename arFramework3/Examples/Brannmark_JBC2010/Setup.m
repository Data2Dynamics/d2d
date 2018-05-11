clear all;
close all;
clc;

%% Compile model
arInit
arLoadModel('model_Mif');
arLoadData('TimeCourse',1);
arLoadData('TwoStepTimeCourse',1);
arLoadData('DoseResponse',1);
arCompileAll;

%% Set pre-equilibration
arSteadyState(1,arFindCondition(ar,'TwoStepTimeCourse'),1);
arSimu(true,true,true);

%% Parameter optimization
arFitLHS(100)
arPlotChi2s

%% Visualization
arSimu(true,true,true);
arPlot;


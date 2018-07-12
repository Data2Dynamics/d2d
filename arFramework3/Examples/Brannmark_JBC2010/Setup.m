% This routine runs the parameter estimation for the model of insulin
% signalling introduced in the publication
%
%   Brännmark, C. and Palmer, R. and Glad, S. T. and Cedersund, G. and 
%   Stralfors, P. (2010). Mass and Information Feedbacks through Receptor
%   Endocytosis Govern Insulin Signaling as Revealed Using a Parameter-free
%   Modeling Framework. J. Biol. Chem., 285(26):20171-20179.
%   doi: 10.1074/jbc.M110.106849

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
arSteadyState(1,arFindCondition(ar,'DoseResponse','insulin_dose_1',0),1:length(ar.model.condition));
arSimu(true,true,true);

%% Set parameters
arLoadPars('bestFit')
arFindInputs;
%% Visualization
arSimu(true,true,true);
arPlot;


% This routine runs the parameter estimation for the model of insulin
% signalling introduced in the publication
%
%   Brännmark, C. and Palmer, R. and Glad, S. T. and Cedersund, G. and 
%   Stralfors, P. (2010). Mass and Information Feedbacks through Receptor
%   Endocytosis Govern Insulin Signaling as Revealed Using a Parameter-free
%   Modeling Framework. J. Biol. Chem., 285(26):20171-20179.
%   doi: 10.1074/jbc.M110.106849

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
arSteadyState(1,arFindCondition(ar,'DoseResponse','insulin_dose_1',0),1:length(ar.model.condition));
arSimu(true,true,true);

%% Set parameters
par_names = {'k1a','k1aBasic','k1b','k1c','k1d','k1e','k1f','k1g','k1r','k21','k22','km2','k3','km3','k_IRP_1Step','k_IRSiP_1Step','k_IRSiP_2Step','k_IRSiP_DosR'};
par_values = [0.38928818520 0.01245274400 0.02002245050 0.36351679280 1580.67826494010 0.00004383400 20.07260350370 286.69946480720 3.63277734420 1.67225033020 0.03638161900 32.59423718910 1.62865902310 0.11310739820 152.96316685360 16760.12031409260 8936.21955740500 13740.43217299910];
par_lb = [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1e1 1e1 1e1 1e1];
par_ub = [5e5  5e5  5e5  5e5  5e5  5e5  5e5  5e5  5e5  5e5  5e5  5e5  5e5  5e5  2e5 2e5 2e5 2e5];
for i = 1:length(par_values)
    arSetPars(par_names{i},log10(par_values(i)),[],[],log10(par_lb(i)),log10(par_ub(i)));
end

%% Parameter optimization
arFitLHS(500)
arPlotChi2s

%% Visualization
arSimu(true,true,true);
arPlot;


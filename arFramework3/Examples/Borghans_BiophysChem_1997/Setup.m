clear all
close all
clc

%% Compile model
arInit
arLoadModel('model_Borghans_BiophysChem_1997');
% This model is the second model proposed in the manuscript by Borghans
% et al. (Biophys Chem, 1997). It is also available from the Biomodels
% database (identifier: BIOMD0000000044).
arLoadData('dataset_Borghans_BiophysChem_1997',1);
% This dataset is a trajectory extracted from Figure 1b (right) of the
% manuscript by Borghans et al. (Biophys Chem, 1997)
arCompileAll;

%% Parameter settings
arSetPars('beta_par',log10(1),1,1,-3,5);
arSetPars('v0',log10(2.0),1,1,-3,5);
arSetPars('v1',log10(1.0),1,1,-3,5);
arSetPars('Vm2',log10(6.5),1,1,-3,5);
arSetPars('K2',log10(0.1),1,1,-3,5);
arSetPars('Vm3',log10(19.5),1,1,-3,5);
arSetPars('Ka',log10(0.2),1,1,-3,5);
arSetPars('Ky',log10(0.2),1,1,-3,5);
arSetPars('Kz',log10(0.3),1,1,-3,5);
arSetPars('Kf',log10(1.0),1,1,-3,5);
arSetPars('K_par',log10(10.0),1,1,-3,5);
arSetPars('Vp',log10(2.5),1,1,-3,5);
arSetPars('Vd',log10(80.0),1,1,-3,5);
arSetPars('Kp',log10(1.0),1,1,-3,5);
arSetPars('Kd',log10(0.4),1,1,-3,5);
arSetPars('n_par',log10(4.4),1,1,-3,5);
arSetPars('epsilon_par',log10(0.1),1,1,-3,5);

arSetPars('scale',log10(1),1,1,-3,5);
arSetPars('offset',log10(0.2),1,1,-3,5);
arSetPars('sigma',log10(0.2),1,1,-3,5);

% Best fit found using a sequence of local optimization problems
% staring with the optimization of the initial conditions and the
% measurement parameters
arLoadPars('20181030T124628_Best');

arFitLHS(1)
arPlotChi2s

arSimu
arPlot

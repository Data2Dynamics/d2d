clear all
close all
clc

%% Compile model
arInit
arLoadModel('model_Elowitz_Nature2000');
arLoadData('dataset_Elowitz_Nature2000',1);
% This dataset is a trajectory extracted from Figure 3b (blue curve) of the
% manuscript by Elowitz & Leibler, Nature, 403, 2000
arCompileAll;

%% Parameter settings
% Parameters of repressilator based on manuscript by
% Elowitz & Leibler, Nature, 403, 2000
arSetPars('n_Hill' ,log10(2),0);
arSetPars('KM',log10(40),1,1,-3,5);
arSetPars('tps_active' ,log10(5e-1),1,1,-3,5);
arSetPars('tps_repr' ,log10(5e-4),1,1,-3,5);
arSetPars('tau_mRNA' ,log10(2.0),1,1,-2,3);
arSetPars('tau_prot',log10(10.0),1,1,-2,3);
arSetPars('eff',log10(20),1,1,0,3);

% Parameters for reporter (selected here)
arSetPars('tau_mRNA_GFP' ,log10(2.0),1,1,-2,3);
arSetPars('tau_prot_GFP',log10(90.0) ,1,1,0,4); % Lutz & Bujard, NAR, 25, 1997.
arSetPars('eff_GFP',log10(1),0); % fixed to remove non-identifiability of eff_GFP*scale
arSetPars('scale',log10(1),1,1,-3,5); 
arSetPars('background',log10(1e-1),1,1,-2,2);

arFitLHS(500)
arPlotChi2s

arPlot

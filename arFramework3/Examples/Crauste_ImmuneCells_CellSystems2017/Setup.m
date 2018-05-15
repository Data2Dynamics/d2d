clear all
close all
clc

% Load models & data
arInit;
arLoadModel('ModelImmuneCells');
arLoadData('DataCrausteImmune');
arCompileAll;

ar.config.fiterrors = -1;

ar.p(5:end) = log10([ ...   
    0.59; % delta_EL 
    0.025; % delta_LM
    0.009; % delta_NE 
    21.5e-6; % mu_EE 
    3.6e-8; % mu_LE 
    7.5e-6; % mu_LL
    0.75; % mu_N 
    5.5e-2; % mu_P 
    1.8e-7; % mu_PE 
    1.8e-5; % mu_PL 
    0.64; % rho_E 
    0.15; % rho_P
    ]); 

ar.lb([9 10 13]) = -10;

%Set fmincon interior-point optimizer
ar.config.maxsteps = 1.e4;
ar.config.optimizer = 2;
ar.config.atol = 1.e-10;
ar.config.rtol = 1.e-10;
ar.config.optim.MaxIter = 1.e4;

arFit
arPlot 


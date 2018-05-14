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
    0.587811356745978;%delta_EL
    0.0249390061548269 % delta_LM
    0.00941887996782282; % delta_NEit 
    2.15367872500831e-05; % mu_EE
    3.57810426702153e-08 % mu_LE
    7.49479597712756e-06; % mu_LL
    0.746532818498732; % mu_N
    0.0551056351722066; % mu_P
    1.79033591694338e-07; % mu_PE
    1.82649004874779e-05; % mu_PL
    0.642885381731727; % % rho_E
    0.150783084474492 % rho_P
    ]);

ar.lb(:) = -10;
ar.config.maxsteps = 1e6;
arFitLHS(100)
arPlot

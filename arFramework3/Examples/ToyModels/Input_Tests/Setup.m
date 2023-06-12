% close all; clear all; clearvars -global;
% Load models & data

arInit;
arLoadModel('input_tests');
arCompileAll(true);

%% set parameters

% injection amount, here 5 nM
arSetPars('bolus_amount', 5, 2, 0, 0, 10);

% injection time point, here at 40 min
arSetPars('injection_timepoint', 40, 2, 0, 0, 50);

% injection duration, here at 1 min
arSetPars('injection_duration', 1, 2, 0, 0, 10);

% step settings
arSetPars('pre_step', 0.5, 0, 0);
arSetPars('post_step', 1.5, 0, 0);

% Size of each step for the double ones
arSetPars('step_size', .5, 0, 0);

% Smoothness of the steps
arSetPars('smoothness', 1.0, 0, 0);
arFindInputs;
arPlot2;

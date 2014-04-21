% Load models & data

arInit;
arLoadModel('bolus_test');
arCompileAll;

%% set parameters for bolus injection

% injection amount, here 5 nM
arSetPars('bolus_amount', 5, 2, 0, 0, 10);

% injection time point, here at 40 min
arSetPars('injection_timepoint', 40, 2, 0, 0, 10);

% injection duration, here at 1 min
arSetPars('injection_duration', 1, 2, 0, 0, 10);

arPlot;

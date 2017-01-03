% Load models & data

arInit;
arLoadModel('wash_bolus_test');
arCompileAll;

%% set parameters for bolus injection

arSetPars('bolus1_amount', 5, 2, 0, 0, 10);             % first injection amount, here 5 nM
arSetPars('injection1_timepoint', 30, 2, 0, 0, 50);     % first injection time point, here at 40 min
arSetPars('injection1_duration', 1, 2, 0, 0, 10);       % first injection duration, here at 1 min

arSetPars('wash_strength', 10, 2, 0, 0, 10);            % washing strength, choose large enough
arSetPars('wash_timepoint', 50, 2, 0, 0, 100);          % washing time point, here at 40 min
arSetPars('wash_duration', 1, 2, 0, 0, 10);             % washing duration, here at 1 min

arSetPars('bolus2_amount', 3, 2, 0, 0, 10);             % second injection amount, here 5 nM
arSetPars('injection2_timepoint', 60, 2, 0, 0, 100);    % second injection time point, here at 40 min
arSetPars('injection2_duration', 1, 2, 0, 0, 10);       % second injection duration, here at 1 min


arPlot;

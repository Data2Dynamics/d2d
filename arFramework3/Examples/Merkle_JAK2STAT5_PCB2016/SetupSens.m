%% Load models & data

arInit;

arLoadModel('jak2_stat5_cfue_sens');
arLoadModel('jak2_stat5_h838_l1_final_sens');

arCompileAll;

arLoadPars('sensitivityanalysis')
ar.lb = 10.^ar.lb;
ar.ub = 10.^ar.ub;

% arSetPars(pLabel, p, qFit, qLog10, lb, ub, type, meanp, stdp)
arSetPars('ActD',0,0,0,0,0)
arSetPars('SOCS3oe',0,0,0,0,0)
arSetPars('overexp',0,0,0,0,0)

arSetPars('epo_level',-6,0,1,-8,-4)

arFindInputs


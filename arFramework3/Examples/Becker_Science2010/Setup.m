% Load models & data

arInit;

arLoadModel('epo_int_rep');
arLoadModel('epo_binding');
arLoadData('Epo_alpha_BaF3_Exp1_cpm_rep','epo_int_rep');
arLoadData('Epo_binding_rep','epo_binding');

arCompileAll;

% parameters with prior information
arSetPars('init_Epo', log10(2100), 1, 1, -5, 4, 1, log10(2100), 0.02);

% load best fit parameter values
arLoadPars('BestFit');

arPlot;
arPrint;
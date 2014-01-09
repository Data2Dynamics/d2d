% Load models & data

arInit;

arLoadModel('epo_int_rep',1);
arLoadModel('epo_binding',2);
arLoadData('Epo_alpha_BaF3_Exp1_cpm_rep',1);
arLoadData('Epo_binding_rep',2);

arCompileAll(true);

% parameters with prior information
arSetPars('init_Epo', log10(2100), 1, 1, -5, 4, 1, log10(2100), 0.02);

% load best fit parameter values
arLoadPars('BestFit');

arPlot;
arPrint;
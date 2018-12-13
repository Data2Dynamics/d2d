% arFitLHSMerge
%
% merge previous fit sequence into current workspace.
% WARNING, fit setting have to be the same. This will not be double checked
% 
% arFitLHS(10)
% arSave('Fit1')
% arFitLHS(90)
% arFitLHSMerge
% -> arStruct with 100 fits

function arFitLHSMerge

global ar

[~, filename] = fileChooser('./Results', true);

fname = ['./Results/' filename '/workspace.mat'];
S = load(fname);

ar.chi2s = [ar.chi2s S.ar.chi2s];
ar.chi2sconstr = [ar.chi2sconstr S.ar.chi2sconstr];
ar.chi2s_start = [ar.chi2s_start S.ar.chi2s_start];
ar.chi2sconstr_start = [ar.chi2sconstr_start S.ar.chi2sconstr_start];

ar.ps = [ar.ps; S.ar.ps];
ar.ps_errors = [ar.ps_errors; S.ar.ps_errors];
ar.ps_start = [ar.ps_start; S.ar.ps_start];

ar.exitflag = [ar.exitflag S.ar.exitflag];
ar.timing = [ar.timing S.ar.timing];
ar.fun_evals = [ar.fun_evals S.ar.fun_evals];
ar.optim_crit = [ar.optim_crit S.ar.optim_crit];
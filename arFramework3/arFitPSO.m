% run Merrimack PSO algorithm

function [pFit, chi2, resnorm, exitflag, output, lambda, jac] = ...
    arFitPSO(LB, UB)

global ar
global fit

% default options
PSOopts = mnb_PSOOptions();

% Objective function options
PSOopts.ObjFuncData = ar;
PSOopts.ObjFuncDir  = '';

PSOopts.FN_OptState = [];
PSOopts.FN_IntState = [];

PSOopts.MaxIter = ar.config.optim.MaxIter;
PSOopts.Verbose = ~strcmp(ar.config.optim.Display, 'off');

% % Becker et al
% PSOopts.NumSwp      = 100;
% PSOopts.lmin        = 4;
% PSOopts.wcount_min  = 3;
% PSOopts.wcount_max  = 6;
% % PSOopts.MaxIter     = 5000;
% % PSOopts.MaxIter     = 100;

% Bachmann et al
PSOopts.NumSwp      = 160;
PSOopts.lmin        = 4;
PSOopts.wcount_min  = 3;
PSOopts.wcount_max  = 5;
% PSOopts.MaxIter     = 20000;

[pFit, chi2, ExitReason, gf_bestHis] = mnb_PSOFit('arFitPSOFkt',LB,UB,PSOopts);

if(strcmp(ExitReason, 'MaxIter Reached'))
    exitflag = 50;
elseif(strcmp(ExitReason, 'FitCount Reache'))
    exitflag = 51;
else
    exitflag = 0;
end

fit.chi2_hist = gf_bestHis';

resnorm = [];
output = struct([]);
output(1).iterations = nan;
lambda = [];
jac = [];

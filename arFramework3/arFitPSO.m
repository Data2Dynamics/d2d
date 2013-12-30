% run Merrimack PSO algorithm

function [pFit, chi2, resnorm, exitflag, output, lambda, jac] = ...
    arFitPSO(LB, UB, ar)

PSOopts = mnb_PSOOptions();

PSOopts.NumSwp      = 100;
PSOopts.lmin        = 4;
PSOopts.wcount_min  = 3;
PSOopts.wcount_max  = 6;


% Objective function options
PSOopts.ObjFuncData = ar;
PSOopts.ObjFuncDir  = '';

PSOopts.MaxIter = ar.config.optim.MaxIter;
PSOopts.Verbose = ~strcmp(ar.config.optim.Display, 'off');

[pFit, chi2, ExitReason] = mnb_PSOFit('arFitPSOFkt',LB,UB,PSOopts);

if(strcmp(ExitReason, 'MaxIter Reached'))
    exitflag = 50;
elseif(strcmp(ExitReason, 'FitCount Reache'))
    exitflag = 51;
else
    exitflag = 0;
end

resnorm = [];
output = struct([]);
output(1).iterations = nan;
lambda = [];
jac = [];

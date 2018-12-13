% arLocalLHS(nfit, [dp], [randomseed])
% 
% Local LHS fitting with nfit inintial guesses in the region ar.p +/- dp.
% 
%   nfit        number of fits
%   dp          bounds/distance of optimal parameter set [0.01]
%   randomseed  which seed, rng(randomseed) []
%
% This funciton might be usefull to check convergence of optimization.
% The original minimum should be found repeatedly.
% 
% Example:
% arFit
% arLocalLHS(100)
% plot(sort(ar.chi2s),'.');
%
% See also arFitLHS

function arLocalLHS(nfit, dp, randomseed)
if(~exist('dp','var') || isempty(dp))
    dp = 0.01;
end
if(~exist('randomseed','var') || isempty(randomseed))
    randomseed = [];
end
global ar

log_fit_history = false;
backup_save = false;


ub = ar.ub + 0.0;
lb = ar.lb + 0.0;

ar.lb = max(ar.lb,ar.p-dp);
ar.ub = min(ar.ub,ar.p+dp);

try
    ps = arRandomPars(nfit, randomseed);
    arFits(ps, log_fit_history, backup_save);
catch ERR
    ar.ub = ub;
    ar.lb = lb;
    rethrow(ERR)
end

ar.ub = ub;
ar.lb = lb;




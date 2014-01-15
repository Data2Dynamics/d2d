% fit sequence using latin hyper cube sampling
% run on MATLAB cluster
%
% arFitLHSCluster(n, randomseed, log_fit_history)
%
% n:                number of runs      [10]
% randomseed:                           rng(randomseed)
% log_fit_history                       [false]

function arFitLHSCluster(n, randomseed, log_fit_history)
if(~exist('n','var'))
    n = 10;
end
if(~exist('randomseed','var'))
    rng('shuffle');
    rngsettings = rng;
    randomseed = rngsettings.Seed;
end
if(~exist('log_fit_history','var'))
    log_fit_history = false;
end
arFitLHS(n, randomseed, log_fit_history, false, true);
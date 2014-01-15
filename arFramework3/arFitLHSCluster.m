% fit sequence using latin hyper cube sampling
% run on MATLAB cluster
%
% arFitLHSCluster(n, randomseed, log_fit_history)
%
% n:                number of runs      [10]
% randomseed:                           rng(randomseed)
% log_fit_history                       [false]

function arFitLHSCluster(n, randomseed, log_fit_history)
arFitLHS(n, randomseed, log_fit_history, backup_save, true);
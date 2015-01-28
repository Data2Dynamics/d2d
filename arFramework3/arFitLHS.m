% fit sequence using 
%   - latin hyper cube sampling (ar.config.useLHS = true)
%   - random sampling from prior
%
% arFitLHS(n, randomseed, log_fit_history, backup_save, use_cluster)
%
% n:                number of runs      [10]
% randomseed:                           rng(randomseed)
% log_fit_history                       [false]
% backup_save                           [false]
% use_cluster                           [false]

function arFitLHS(n, randomseed, log_fit_history, backup_save, use_cluster)

if(~exist('n','var'))
    n = 10;
end
if(~exist('randomseed','var'))
    randomseed = [];
end
if(~exist('log_fit_history','var'))
    log_fit_history = false;
end
if(~exist('backup_save','var'))
    backup_save = false;
end
if(~exist('use_cluster','var'))
    use_cluster = false;
end

% generate random values
ps = arRandomPars(n, randomseed);

if(~use_cluster)
    arFits(ps, log_fit_history, backup_save);
else
    arFitsCluster(ps, log_fit_history);
end



% fit sequence using latin hyper cube sampling
%
% arFitLHS(n, randomseed, log_fit_history, backup_save, use_cluster)
%
% n:                number of runs      [10]
% randomseed:                           rng(randomseed)
% log_fit_history                       [false]
% backup_save                           [false]
% use_cluster                           [false]

function arFitLHS(n, randomseed, log_fit_history, backup_save, use_cluster)

global ar

if(~exist('n','var'))
    n = 10;
end
if(exist('rng','file')~=0)
    if(exist('randomseed','var'))
        ar.lhs_seed = randomseed;
        rng(randomseed);
    else
        rng('shuffle');
        rngsettings = rng;
        ar.lhs_seed = rngsettings.Seed;
    end
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

ps = ones(n,1) * ar.p;

q_select = ar.qFit==1;
psrand = lhsdesign(n,sum(q_select));
psrand = psrand .* (ones(n,1)*(ar.ub(q_select) - ar.lb(q_select)));
psrand = psrand + (ones(n,1)*ar.lb(q_select));

ps(:,q_select) = psrand;

if(~use_cluster)
    arFits(ps, log_fit_history, backup_save);
else
    arFitsCluster(ps, log_fit_history);
end



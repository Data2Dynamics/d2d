% fit sequence using latin hyper cube sampling
%
% arFitLHS(n, append, randomseed, dynamic_only, plot_summary, log_fit_history, backup_save)
%
% n:        number of runs      [10]
% append:                       [false]
% randomseed:                   rng(randomseed)
% dynamic_only                  [false]
% plot_summary                  [false]
% log_fit_history               [false]
% backup_save                   [false]

function arFitLHS(n, append, randomseed, dynamic_only, plot_summary, log_fit_history, backup_save)

global ar

if(~exist('n','var'))
    n = 10;
end
if(~exist('append','var'))
    append = false;
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
if(~exist('dynamic_only','var'))
    dynamic_only = false;
end
if(~exist('plot_summary','var'))
    plot_summary = false;
end
if(~exist('log_fit_history','var'))
    log_fit_history = false;
end
if(~exist('backup_save','var'))
    backup_save = false;
end

ps = ones(n,1) * ar.p;

if(dynamic_only)
    q_select = ar.qFit==1 & ar.qDynamic==1;
else
    q_select = ar.qFit==1;
end

psrand = lhsdesign(n,sum(q_select));
psrand = psrand .* (ones(n,1)*(ar.ub(q_select) - ar.lb(q_select)));
psrand = psrand + (ones(n,1)*ar.lb(q_select));

ps(:,q_select) = psrand;

arFits(ps, append, dynamic_only, plot_summary, log_fit_history, backup_save);



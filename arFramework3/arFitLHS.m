% fit sequence using latin hyper cube sampling
%
% arFitLHS(n, randomseed, log_fit_history)
%
% n:                number of runs      [10]
% randomseed:                           rng(randomseed)
% log_fit_history                       [false]

function arFitLHS(n, randomseed, log_fit_history)

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

ps = ones(n,1) * ar.p;

q_select = ar.qFit==1;
psrand = lhsdesign(n,sum(q_select));
psrand = psrand .* (ones(n,1)*(ar.ub(q_select) - ar.lb(q_select)));
psrand = psrand + (ones(n,1)*ar.lb(q_select));

ps(:,q_select) = psrand;

arFits(ps, log_fit_history);



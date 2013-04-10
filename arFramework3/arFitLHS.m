% fit sequence using latin hyper cube sampling
%
% arFitLHS(n, append, randomseed, dynamic_only)
%
% n:        number of runs      [10]
% append:                       [false]
% randomseed:   rng(randomseed)
% dynamic_only                  [true]

function arFitLHS(n, append, randomseed, dynamic_only, plot_summary)

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
        ar.lhs_seed = now;
        rng(ar.lhs_seed);
    end
end
if(~exist('dynamic_only','var'))
    dynamic_only = true;
end
if(~exist('plot_summary','var'))
    plot_summary = false;
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

arFits(ps, append, dynamic_only, plot_summary);

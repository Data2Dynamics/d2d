% chi2 sequence using latin hyper cube sampling
%
% arChi2LHS(n, sensis, randomseed, silent, use_cluster)
%
% n:            number of runs          [10]
% sensis:       use sensitivities       [false]
% randomseed:                           rng(randomseed)
% silent:       no output               [false]
% use_cluster:                          [false]

function arChi2LHS(n, sensis, randomseed, silent, use_cluster)

global ar

if(~exist('n','var'))
    n = 10;
end
if(~exist('sensis','var'))
    sensis = false;
end
if(exist('rng','file')~=0)
    if(exist('randomseed','var') && ~isempty(randomseed))
        ar.lhs_seed = randomseed;
        rng(randomseed);
    else
        rng('shuffle');
        rngsettings = rng;
        ar.lhs_seed = rngsettings.Seed;
    end
end
if(~exist('silent','var'))
    silent = false;
end
if(~exist('use_cluster','var'))
    use_cluster = false;
end

ps = ones(n,1) * ar.p;
psrand = lhsdesign(n,sum(ar.qFit==1));

psrand = psrand .* (ones(n,1)*(ar.ub(ar.qFit==1) - ar.lb(ar.qFit==1)));
psrand = psrand + (ones(n,1)*ar.lb(ar.qFit==1));

ps(:,ar.qFit==1) = psrand;

if(~use_cluster)
    arChi2s(ps, sensis, silent);
else
    arChi2sCluster(ps, sensis, silent);
end


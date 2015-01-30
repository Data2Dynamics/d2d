% chi2 sequence using latin hyper cube sampling
%   - latin hyper cube sampling (ar.config.useLHS = true)
%   - random sampling from prior
%
% arChi2LHS(n, sensis, randomseed, silent, use_cluster)
%
% n:            number of runs          [10]
% sensis:       use sensitivities       [false]
% randomseed:                           rng(randomseed)
% silent:       no output               [false]
% use_cluster:                          [false]

function arChi2LHS(n, sensis, randomseed, silent, use_cluster)

if(~exist('n','var'))
    n = 10;
end
if(~exist('sensis','var'))
    sensis = false;
end
if(~exist('randomseed','var'))
    randomseed = [];
end
if(~exist('silent','var'))
    silent = false;
end
if(~exist('use_cluster','var'))
    use_cluster = false;
end

% generate random values
ps = arRandomPars(n, randomseed);

if(~use_cluster)
    arChi2s(ps, sensis, silent);
else
    arChi2sCluster(ps, sensis, silent);
end


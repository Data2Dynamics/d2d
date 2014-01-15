% chi2 sequence using latin hyper cube sampling
% run on MATLAB cluster
%
% arChi2LHSCluster(n, sensis, silent, use_cluster)
%
% n:            number of runs          [10]
% sensis:       use sensitivities       [false]
% silent:       no output               [false]

function arChi2LHSCluster(n, sensis, silent)
if(~exist('n','var'))
    n = 10;
end
if(~exist('sensis','var'))
    sensis = false;
end
if(~exist('silent','var'))
    silent = false;
end
arChi2LHS(n, sensis, silent, true);
% chi2 sequence using 
%   - latin hyper cube sampling (ar.config.useLHS = true)
%   - random sampling from prior
% run on MATLAB cluster
%
% job = arChi2LHSCluster(cluster, clusterpath, pool_size, n, sensis, randomseed, silent)
%
% cluster:      MATLAB cluster object       (see help parcluster)
% clusterpath:  execution path on cluster   ['.']
% pool_size:    additional workers          [ceil(length(cluster.IdleWorkers)/2)]
% n:            number of runs              [10]
% sensis:       use sensitivities           [false]
% randomseed:                               rng(randomseed)
% silent:       show output                 [false]

function varargout = arChi2LHSCluster(cluster, clusterpath, pool_size, n, sensis, randomseed, silent)

global ar
global ar_chi2lhs_cluster

if(isempty(ar_chi2lhs_cluster)) % new job
    if(~exist('clusterpath','var') || isempty(clusterpath))
        clusterpath = '.';
    end
    if(~exist('pool_size','var') || isempty(pool_size))
        pool_size = ceil(length(cluster.IdleWorkers)/2);
    end
    
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

    fprintf('arChi2LHSCluster sending job...');
    ar_chi2lhs_cluster = batch(cluster, @arChi2LHSClusterFun, 1, ...
        {ar, n, sensis, randomseed, silent}, ...
        'CaptureDiary', true, ...
        'CurrentFolder', clusterpath, ...
        'Matlabpool', pool_size);
    fprintf('done\n');
    
elseif(isa(ar_chi2lhs_cluster,'parallel.job.MJSIndependentJob') || ...
        isa(ar_chi2lhs_cluster,'parallel.job.MJSCommunicatingJob') || ...
        isa(ar_chi2lhs_cluster,'parallel.job.CJSIndependentJob') || ...
        isa(ar_chi2lhs_cluster,'parallel.job.CJSCommunicatingJob')) % old job
    if(nargout>0)
        varargout{1} = ar_chi2lhs_cluster;
        return
    end
    try
        fprintf('arChi2LHSCluster (ID %i) ', ar_chi2lhs_cluster.ID);
    catch
        fprintf('arChi2LHSCluster invalid job ID, deleting...\n');
        clear global ar_chi2lhs_cluster
        return
    end
    
    if(~strcmp(ar_chi2lhs_cluster.State, 'finished')) % still running
        fprintf('status %s...\n', ar_chi2lhs_cluster.State);
    else % finished
        try
            fprintf('retrieving results %s...\n', ar_chi2lhs_cluster.State);
            diary(ar_chi2lhs_cluster);
            S = fetchOutputs(ar_chi2lhs_cluster);
            ar = S{1};
            delete(ar_chi2lhs_cluster);
            clear global ar_chi2lhs_cluster
        catch err_id
            delete(ar_chi2lhs_cluster);
            clear global ar_chi2lhs_cluster
            rethrow(err_id);
        end
    end
else
    error('arChi2LHSCluster global variable ar_chi2lhs_cluster is invalid!\n');
end

function ar = arChi2LHSClusterFun(ar2, n, sensis, randomseed, silent)
global ar %#ok<REDEF>
ar = ar2; %#ok<NASGU>
arChi2LHS(n, sensis, randomseed, silent, true);

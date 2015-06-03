% fit sequence using 
%   - latin hyper cube sampling (ar.config.useLHS = true)
%   - random sampling from prior
% run on MATLAB cluster
%
% job = arFitLHSCluster(cluster, clusterpath, pool_size, n, randomseed, log_fit_history)
%
% cluster:          MATLAB cluster object       (see help parcluster)
% clusterpath:      execution path on cluster   ['.']
% pool_size:        additional workers          [ceil(length(cluster.IdleWorkers)/2)]
% n:                number of runs              [10]
% randomseed:                                   rng(randomseed)
% log_fit_history                               [false]

function varargout = arFitLHSCluster(cluster, clusterpath, pool_size, n, randomseed, log_fit_history)

global ar
global ar_fitlhs_cluster

if(isempty(ar_fitlhs_cluster)) % new job
    if(nargin==0)
        error('specify cluster!');
    end
    if(nargout>0)
        error('no job available');
    end
    
    if(~exist('clusterpath','var') || isempty(clusterpath))
        clusterpath = '.';
    end
    if(~exist('pool_size','var') || isempty(pool_size))
        pool_size = ceil(length(cluster.IdleWorkers)/2);
    end
    
    if(~exist('n','var'))
        n = 10;
    end
    if(~exist('randomseed','var'))
        randomseed = [];
    end
    if(~exist('log_fit_history','var'))
        log_fit_history = false;
    end

    fprintf('arFitLHSCluster sending job...');
    ar_fitlhs_cluster = batch(cluster, @arFitLHSClusterFun, 1, ...
        {ar, n, randomseed, log_fit_history}, ...
        'CaptureDiary', true, ...
        'CurrentFolder', clusterpath, ...
        'pool', pool_size, ...
        'AttachedFiles', {'arFitPSOFkt'});
    
    fprintf('done\n');
    
elseif(isa(ar_fitlhs_cluster,'parallel.job.MJSIndependentJob') || ...
        isa(ar_fitlhs_cluster,'parallel.job.MJSCommunicatingJob') || ...
        isa(ar_fitlhs_cluster,'parallel.job.CJSIndependentJob') || ...
        isa(ar_fitlhs_cluster,'parallel.job.CJSCommunicatingJob')) % old job
    if(nargout>0)
        varargout{1} = ar_fitlhs_cluster;
        return
    end
    try
        fprintf('arFitLHSCluster (ID %i) ', ar_fitlhs_cluster.ID);
    catch
        fprintf('arFitLHSCluster invalid job ID, deleting...\n');
        clear global ar_fitlhs_cluster
        return
    end
    
    if(~strcmp(ar_fitlhs_cluster.State, 'finished')) % still running
        fprintf('status %s...\n', ar_fitlhs_cluster.State);
    else % finished
        fprintf('retrieving results %s...\n', ar_fitlhs_cluster.State);
        diary(ar_fitlhs_cluster);
        S = fetchOutputs(ar_fitlhs_cluster);
        ar = S{1};
        delete(ar_fitlhs_cluster);
        clear global ar_fitlhs_cluster
    end
else
    error('arFitLHSCluster global variable ar_fitlhs_cluster is invalid!\n');
end

function ar = arFitLHSClusterFun(ar2, n, randomseed, log_fit_history)
global ar %#ok<REDEF>
ar = ar2; %#ok<NASGU>
arFitLHS(n, randomseed, log_fit_history, false, true);

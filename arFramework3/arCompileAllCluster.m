% compile on a matlab cluster worker
%
% job = arCompileAllCluster(cluster, clusterpath, pool_size)
%
% cluster:      MATLAB cluster object       (see help parcluster)
% clusterpath:  execution path on cluster   ['.']
% pool_size:    additional workers          [ceil(length(cluster.IdleWorkers)/2)]

function varargout = arCompileAllCluster(cluster, clusterpath, pool_size)

global ar
global ar_compileall_cluster

if(isempty(ar_compileall_cluster) || strcmp(ar_compileall_cluster.Status, 'deleted')) % new job
    fprintf('arCompileAllCluster sending job...');
   
    if(~exist('pool_size','var'))
        pool_size = ceil(length(cluster.IdleWorkers)/2);
    end
    
    ar_compileall_cluster = batch(cluster, @arCompileAllClusterFun, 1, ...
        {ar}, ...
        'CaptureDiary', true, ...
        'CurrentFolder', clusterpath, ...
        'pool', pool_size);
    fprintf('done\n');
    
elseif(isa(ar_compileall_cluster,'parallel.job.MJSIndependentJob') || ...
        isa(ar_compileall_cluster,'parallel.job.MJSCommunicatingJob') || ...
        isa(ar_compileall_cluster,'parallel.job.CJSIndependentJob') || ...
        isa(ar_compileall_cluster,'parallel.job.CJSCommunicatingJob')) % old job
    if(nargout>0)
        varargout{1} = ar_compileall_cluster;
        return
    end
    try
        fprintf('arCompileAllCluster (ID %i) ', ar_compileall_cluster.ID);
    catch
        fprintf('arCompileAllCluster invalid job ID, deleting...\n');
        clear global ar_compileall_cluster
        return
    end
    
    if(~strcmp(ar_compileall_cluster.State, 'finished')) % still running
        fprintf('status %s...\n', ar_compileall_cluster.State);
    else % finished
        try
            fprintf('retrieving results %s...\n', ar_compileall_cluster.State);
            diary(ar_compileall_cluster);
            S = fetchOutputs(ar_compileall_cluster);
            ar = S{1};
            delete(ar_compileall_cluster);
            clear global ar_compileall_cluster
        catch err_id
            delete(ar_compileall_cluster);
            clear global ar_compileall_cluster
            rethrow(err_id);
        end
    end
else
    error('arCompileAllCluster global variable ar_compileall_cluster is invalid!\n');
end

function ar = arCompileAllClusterFun(ar2)
global ar %#ok<REDEF>
ar = ar2; %#ok<NASGU>
addpath('/data/Andreas/arFramework3');
arCompileAll;

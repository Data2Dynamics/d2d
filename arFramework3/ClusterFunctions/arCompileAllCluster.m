% job = arCompileAllCluster(cluster, clusterpath, pool_size)
% compile on a matlab cluster worker
%
% cluster:      MATLAB cluster object       (see help parcluster)
% clusterpath:  execution path on cluster   ['.']
% pool_size:    additional workers          [ceil(length(cluster.IdleWorkers)/2) or (cluster.NumWorkers-1)]

function varargout = arCompileAllCluster(cluster, clusterpath, pool_size, addpath)

global ar
global ar_compileall_cluster

if(isempty(ar_compileall_cluster) || strcmp(ar_compileall_cluster.State, 'deleted')) % new job
    fprintf('arCompileAllCluster sending job...');
   
    if(~exist('clusterpath','var') || isempty(clusterpath))
        clusterpath = '.';
    end
    if(~exist('addpath','var') || isempty(addpath))
        addpath = '.';
    end
    if(~exist('pool_size','var') || isempty(pool_size))
        if isfield(cluster,'IdleWorkers') % only exists for certain cluster objects
            pool_size = ceil(length(cluster.IdleWorkers)/2);
        else
            pool_size = cluster.NumWorkers-1;
        end
    end
    
    ar_compileall_cluster = batch(cluster, @arCompileAllClusterFun, 1, ...
        {ar}, ...
        'CaptureDiary', true, ...
        'CurrentFolder', clusterpath, ...
        'AdditionalPaths', {addpath [addpath '/Ccode']}, ...
        'pool', pool_size);
    fprintf('done\n');
    if(nargout>0)
        varargout{1} = ar_compileall_cluster;
        return
    end
    
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
addpath(ar.info.ar_path); %#ok<NODEF>
arCompileAll;

% fit sequence using latin hyper cube sampling
% run on MATLAB cluster
%
% job = arFitLHSCluster(cluster, n, randomseed, log_fit_history)
%
% cluster:          MATLAB cluster object   (see help parcluster)
% n:                number of runs          [10]
% randomseed:                               rng(randomseed)
% log_fit_history                           [false]

function varargout = arFitLHSCluster(cluster, n, randomseed, log_fit_history, pool_size)

global ar
global ar_fitlhs_cluster

if(isempty(ar_fitlhs_cluster)) % new job
    if(nargin==0)
        error('specify cluster!');
    end
    if(nargout>0)
        error('no job available');
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
    
    if(~exist('pool_size','var'))
        pool_size = ceil(length(cluster.IdleWorkers)/2);
    end

    fprintf('arFitLHSCluster sending job...');
    ar_fitlhs_cluster = batch(cluster, @arFitLHSClusterFun, 1, ...
        {ar, n, randomseed, log_fit_history}, ...
        'CaptureDiary', true, ...
        'CurrentFolder', '.', ...
        'AttachedFiles', {'arFitPSOFkt'}, ...
        'Matlabpool', pool_size);
    
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
        try
            fprintf('retrieving results %s...\n', ar_fitlhs_cluster.State);
            diary(ar_fitlhs_cluster);
            S = fetchOutputs(ar_fitlhs_cluster);
            ar = S{1};
            delete(ar_fitlhs_cluster);
            clear global ar_fitlhs_cluster
        catch err_id
            delete(ar_fitlhs_cluster);
            clear global ar_fitlhs_cluster
            rethrow(err_id);
        end
    end
else
    error('arFitLHSCluster global variable ar_fitlhs_cluster is invalid!\n');
end

function ar = arFitLHSClusterFun(ar2, n, randomseed, log_fit_history)
global ar %#ok<REDEF>
ar = ar2; %#ok<NASGU>
arFitLHS(n, randomseed, log_fit_history, false, true);

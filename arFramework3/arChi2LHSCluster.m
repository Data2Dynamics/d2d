% chi2 sequence using latin hyper cube sampling
% run on MATLAB cluster
%
% job = arChi2LHSCluster(cluster, n, sensis, silent)
%
% cluster:      MATLAB cluster object   (see help parcluster)
% n:            number of runs          [10]
% sensis:       use sensitivities       [false]
% silent:       show output             [false]

function varargout = arChi2LHSCluster(cluster, n, sensis, silent)

global ar
global ar_chi2lhs_cluster

if(isempty(ar_chi2lhs_cluster)) % new job    
    if(~exist('n','var'))
        n = 10;
    end
    if(~exist('sensis','var'))
        sensis = false;
    end
    if(~exist('silent','var'))
        silent = false;
    end

    fprintf('arChi2LHSCluster sending job...');
    ar_chi2lhs_cluster = batch(cluster, @arChi2LHSClusterFun, 1, ...
        {ar, n, sensis, silent}, ...
        'CaptureDiary', true, ...
        'CurrentFolder', '.');
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

function ar = arChi2LHSClusterFun(ar2, n, sensis, silent)
global ar %#ok<REDEF>
ar = ar2; %#ok<NASGU>
arChi2LHS(n, sensis, silent, true);

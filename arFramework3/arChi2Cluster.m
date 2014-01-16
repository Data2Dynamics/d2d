% send a single evaluation to a matlab cluster worker
%
% job = arChi2Cluster(cluster, sensi)
%
% cluster:      MATLAB cluster object       (see help parcluster)
% sensi:        propagate sensitivities     [false]

function varargout = arChi2Cluster(cluster, sensi)

global ar
global ar_chi2_cluster

if(isempty(ar_chi2_cluster)) % new job
    if(~exist('n','var'))
        sensi = false;
    end

    fprintf('arChi2Cluster sending job...');
    ar_chi2_cluster = batch(cluster, @arChi2ClusterFun, 1, {ar, sensi}, ...
        'CaptureDiary', true, ...
        'CurrentFolder', '.');
    fprintf('done\n');
    
elseif(isa(ar_chi2_cluster,'parallel.job.MJSIndependentJob')) % old job
    try
        fprintf('arChi2Cluster (ID %i) ', ar_chi2_cluster.ID);
    catch
        fprintf('arChi2Cluster invalid job ID, deleting...\n');
        clear global ar_chi2_cluster
        return
    end
    
    if(~strcmp(ar_chi2_cluster.State, 'finished')) % still running
        fprintf('status %s...\n', ar_chi2_cluster.State);
    else % finished
        try
            fprintf('retrieving results %s...\n', ar_chi2_cluster.State);
            diary(ar_chi2_cluster);
            S = fetchOutputs(ar_chi2_cluster);
            ar = S{1};
            delete(ar_chi2_cluster);
            clear global ar_chi2_cluster
        catch err_id
            delete(ar_chi2_cluster);
            clear global ar_chi2_cluster
            rethrow(err_id);
        end
    end
else
    error('arChi2Cluster global variable ar_chi2_cluster is invalid!\n');
end

if(nargout>0)
    varargout{1} = ar_chi2_cluster;
end

function ar = arChi2ClusterFun(ar, sensi)
ar = arChi2(ar, sensi);

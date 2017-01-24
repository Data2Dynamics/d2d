% send a single fit to a matlab cluster worker
%
% job = arFitCluster(cluster)
%
% cluster:      MATLAB cluster object   (see help parcluster)

function varargout = arFitCluster(cluster)

global ar
global ar_fit_cluster

if(isempty(ar_fit_cluster)) % new job
    fprintf('arFitCluster sending job...');
    ar_fit_cluster = batch(cluster, @arFitClusterFun, 1, {ar}, ...
        'CaptureDiary', true, ...
        'CurrentFolder', '.');
    fprintf('done\n');
    
elseif(isa(ar_fit_cluster,'parallel.job.MJSIndependentJob') || ...
        isa(ar_fit_cluster,'parallel.job.MJSCommunicatingJob') || ...
        isa(ar_fit_cluster,'parallel.job.CJSIndependentJob') || ...
        isa(ar_fit_cluster,'parallel.job.CJSCommunicatingJob')) % old job
    if(nargout>0)
        varargout{1} = ar_fit_cluster;
        return
    end
    try
        fprintf('arFitCluster (ID %i) ', ar_fit_cluster.ID);
    catch
        fprintf('arFitCluster invalid job ID, deleting...\n');
        clear global ar_fit_cluster
        return
    end
    
    if(~strcmp(ar_fit_cluster.State, 'finished')) % still running
        fprintf('status %s...\n', ar_fit_cluster.State);
    else % finished
        try
            fprintf('retrieving results %s...\n', ar_fit_cluster.State);
            diary(ar_fit_cluster);
            S = fetchOutputs(ar_fit_cluster);
            ar = S{1};
            t_onCluster = datenum(ar_fitlhs_cluster.FinishTime([5:20 25:end]),'mmm dd HH:MM:SS YYYY')...
                    - datenum(ar_fitlhs_cluster.StartTime([5:20 25:end]),'mmm dd HH:MM:SS YYYY');
                fprintf('total time on cluster %ih %im %is \n',str2num(datestr(t_onCluster,'HH;MM;SS')));
            delete(ar_fit_cluster);
            clear global ar_fit_cluster
        catch err_id
            delete(ar_fit_cluster);
            clear global ar_fit_cluster
            rethrow(err_id);
        end
    end
else
    error('arFitCluster global variable ar_fit_cluster is invalid!\n');
end

function ar = arFitClusterFun(ar)
ar = arFit(ar);

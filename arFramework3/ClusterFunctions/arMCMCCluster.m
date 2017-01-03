% send a single mcmc run to a matlab cluster worker
%
% job = arMCMCCluster(cluster, nruns, nburnin, method, append, nthinning)
%
% cluster:      MATLAB cluster object   (see help parcluster)

function varargout = arMCMCCluster(cluster, nruns, nburnin, method, append, nthinning)

global ar
global ar_mcmc_cluster

if(isempty(ar_mcmc_cluster)) % new job
    fprintf('arMCMCCluster sending job...');
    
    if(~exist('nruns','var'))
        nruns = 1000;
    end
    nwindow = sum(ar.qFit == 1)*50;
    if(~exist('nburnin','var') || nburnin == 0)
        nburnin = 0;
        if(method==4)
            nburnin = nwindow * 50;
        end
    end
    if(~exist('method','var'))
        method = 1;
    end
    if(~exist('append','var'))
        append = false;
    end
    if(~exist('nthinning','var'))
        nthinning = 1;
    end
    
    ar_mcmc_cluster = batch(cluster, @arMCMCClusterFun, 1, ...
        {ar, nruns, nburnin, method, append, nthinning}, ...
        'CaptureDiary', true, ...
        'CurrentFolder', '.');
    fprintf('done\n');
    
elseif(isa(ar_mcmc_cluster,'parallel.job.MJSIndependentJob') || ...
        isa(ar_mcmc_cluster,'parallel.job.MJSCommunicatingJob') || ...
        isa(ar_mcmc_cluster,'parallel.job.CJSIndependentJob') || ...
        isa(ar_mcmc_cluster,'parallel.job.CJSCommunicatingJob')) % old job
    if(nargout>0)
        varargout{1} = ar_mcmc_cluster;
        return
    end
    try
        fprintf('arMCMCCluster (ID %i) ', ar_mcmc_cluster.ID);
    catch
        fprintf('arMCMCCluster invalid job ID, deleting...\n');
        clear global ar_mcmc_cluster
        return
    end
    
    if(~strcmp(ar_mcmc_cluster.State, 'finished')) % still running
        fprintf('status %s...\n', ar_mcmc_cluster.State);
    else % finished
        try
            fprintf('retrieving results %s...\n', ar_mcmc_cluster.State);
            diary(ar_mcmc_cluster);
            S = fetchOutputs(ar_mcmc_cluster);
            ar = S{1};
            delete(ar_mcmc_cluster);
            clear global ar_mcmc_cluster
        catch err_id
            delete(ar_mcmc_cluster);
            clear global ar_mcmc_cluster
            rethrow(err_id);
        end
    end
else
    error('arMCMCCluster global variable ar_mcmc_cluster is invalid!\n');
end

function ar = arMCMCClusterFun(ar2, nruns, nburnin, method, append, nthinning)
global ar %#ok<REDEF>
ar = ar2; %#ok<NASGU>
global arWaitbarGlobal
arWaitbarGlobal.showWindow = 0;
arMCMC(nruns, nburnin, method, append, nthinning);

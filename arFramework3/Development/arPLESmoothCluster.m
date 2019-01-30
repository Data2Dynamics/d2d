% smooth profile likelilhood calculation on cluster
%
% job = arPLESmoothCluster(cluster, [clusterpath, pool_size, jk])
%
% cluster:          MATLAB cluster object       (see help parcluster)
% clusterpath:      execution path on cluster   ['.']
% pool_size:        additional workers          [ceil(length(cluster.IdleWorkers)/2) or (cluster.NumWorkers-1)]
% jk:               parameter index or indices  [all fit parameters]

function varargout = arPLESmoothCluster(cluster, clusterpath, pool_size, jk)

global ar
global ar_plecalc_cluster

if(isempty(ar_plecalc_cluster)) % new job
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
        if isfield(cluster,'IdleWorkers') % only exists for certain cluster objects
            pool_size = ceil(length(cluster.IdleWorkers)/2);
        else
            pool_size = cluster.NumWorkers-1;
        end
    end
    
    if(~exist('jk','var'))
        jk = find(ar.qFit==1);
        jk = jk(jk<=length(ar.ple.chi2s));
    end

    fprintf('arPLESmoothCluster sending job...');
    ar_plecalc_cluster = batch(cluster, @arPLESmoothClusterFun, 1, ...
        {ar, jk}, ...
        'CaptureDiary', true, ...
        'CurrentFolder', clusterpath, ...
        'pool', pool_size);
    
    fprintf('done\n');
    
elseif(isa(ar_plecalc_cluster,'parallel.job.MJSIndependentJob') || ...
        isa(ar_plecalc_cluster,'parallel.job.MJSCommunicatingJob') || ...
        isa(ar_plecalc_cluster,'parallel.job.CJSIndependentJob') || ...
        isa(ar_plecalc_cluster,'parallel.job.CJSCommunicatingJob')) % old job
    if(nargout>0)
        varargout{1} = ar_plecalc_cluster;
        return
    end
    try
        fprintf('arPLESmoothCluster (ID %i) ', ar_plecalc_cluster.ID);
    catch
        fprintf('arPLESmoothCluster invalid job ID, deleting...\n');
        clear global ar_plecalc_cluster
        return
    end
    
    if(~strcmp(ar_plecalc_cluster.State, 'finished')) % still running
        fprintf('status %s...\n', ar_plecalc_cluster.State);
    else % finished
        fprintf('retrieving results %s...\n', ar_plecalc_cluster.State);
        diary(ar_plecalc_cluster);
        S = fetchOutputs(ar_plecalc_cluster);
        ar = S{1};
        delete(ar_plecalc_cluster);
        clear global ar_plecalc_cluster
    end
else
    error('arPLESmoothCluster global variable ar_plecalc_cluster is invalid!\n');
end



function ar = arPLESmoothClusterFun(ar2, jk)

global ar %#ok<REDEF>
ar = ar2;  %#ok<*NASGU>

njk = length(jk);

chi2s = cell(1,njk);
errors = cell(1,njk);
ps = cell(1,njk);
run = zeros(1,njk);
chi2Reset = nan(1,njk);

tic;
startTime = clock;
arShowProgressParFor(njk);
ar1 = ar; %#ok<NODEF>
parfor j=1:njk
    ar2 = ar1;
    
    ar2 = arPLESmooth(ar2, jk(j));
    
    chi2s{j} = ar2.ple.chi2s{jk(j)};
    errors{j} = ar2.ple.errors{jk(j)};
    ps{j} = ar2.ple.ps{jk(j)};
    run(j) = ar2.ple.run(jk(j));
    chi2Reset(j) = ar2.ple.chi2Reset(jk(j));

    arShowProgressParFor(j, njk, startTime)
end
arShowProgressParFor(0);
toc;

ar.ple.chi2s(jk) = chi2s;  %#ok<STRNU>
ar.ple.errors(jk) = errors;  %#ok<STRNU>
ar.ple.ps(jk) = ps;   %#ok<STRNU>
ar.ple.run(jk) = run;  %#ok<STRNU>
ar.ple.chi2Reset(jk) = chi2Reset;  %#ok<STRNU>

%% backup save on cluster before transfer
% Sometime the automatic data transfer from the cluster does not work.
% So it is better to save a local copy of the final results that 
% can be retrieved by file transfer.

save('arPLESmoothCluster_backup.mat', 'ar');



% profile likelilhood calculation on cluster
%
% job = arPLECalcCluster(cluster, clusterpath, pool_size, jk, n)
%
% cluster:          MATLAB cluster object       (see help parcluster)
% clusterpath:      execution path on cluster   ['.']
% pool_size:        additional workers          [ceil(length(cluster.IdleWorkers)/2)]
% jk:               parameter index or indices          [all fit parameters]
% n:                number of ple steps up and down     [50]

function varargout = arPLECalcCluster(cluster, clusterpath, pool_size, jk, n)

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
        pool_size = ceil(length(cluster.IdleWorkers)/2);
    end
    
    if(~exist('jk','var'))
        jk = find(ar.qFit==1);
    end
    if(~exist('n','var'))
        n = 50;
    end

    fprintf('arPLECalcCluster sending job...');
    ar_plecalc_cluster = batch(cluster, @arPLECalcClusterFun, 1, ...
        {ar, jk, n}, ...
        'CaptureDiary', true, ...
        'CurrentFolder', clusterpath, ...
        'pool', pool_size, ...
        'AttachedFiles', {'arFitPSOFkt'});
    
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
        fprintf('arPLECalcCluster (ID %i) ', ar_plecalc_cluster.ID);
    catch
        fprintf('arPLECalcCluster invalid job ID, deleting...\n');
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
    error('arPLECalcCluster global variable ar_plecalc_cluster is invalid!\n');
end



function ar = arPLECalcClusterFun(ar2, jk, n)

global ar %#ok<REDEF>
ar = ar2;  %#ok<*NASGU>

if(~isfield(ar, 'ple')) %#ok<*NODEF>
    ar.ple.chi2s = {}; %#ok<*STRNU>
    ar.ple.errors = {};
    ar.ple.ps = {};
    ar.ple.run = zeros(size(ar.p));
    ar.ple.chi2Reset = nan(size(ar.p));
    ar.ple.alpha = 0.05;
    ar.ple.ndof = 1;
end
ar.ple.pStart = ar.p;

njk = length(jk);

chi2s = cell(1,njk);
errors = cell(1,njk);
ps = cell(1,njk);
run = zeros(1,njk);
chi2Reset = nan(1,njk);

startTime = clock;
arShowProgressParFor(njk);
ar1 = ar;
parfor j=1:njk
    ar2 = ar1;
    
    ar2 = arPLECalc(ar2, jk(j), n);
    
    chi2s{j} = ar2.ple.chi2s{jk(j)};
    errors{j} = ar2.ple.errors{jk(j)};
    ps{j} = ar2.ple.ps{jk(j)};
    run(j) = ar2.ple.run(jk(j));
    chi2Reset(j) = ar2.ple.chi2Reset(jk(j));

    arShowProgressParFor(j, njk, startTime)
end
arShowProgressParFor(0);

ar.ple.chi2s(jk) = chi2s;
ar.ple.errors(jk) = errors;
ar.ple.ps(jk) = ps;
ar.ple.run(jk) = run;
ar.ple.chi2Reset(jk) = chi2Reset;

%% backup save on cluster before transfer
% Sometime the automatic data transfer from the cluster does not work.
% So it is better to save a local copy of the final results that 
% can be retrieved by file transfer.

save('arPLECalcCluster_backup.mat', 'ar');



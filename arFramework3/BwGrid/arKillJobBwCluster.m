function arKillJobBwCluster(jobId, partition)

% arKillJobBwCluster(jobId, [partition])
%
% Kill jobs on slurm queueing system.
%   jobId       Vector of jobIds OR 'all' to kill all of your jobs
%   partition   Only kill jobs on specified partition. Meaningful only when
%               using 'all' for jobid
%
% See also: arHelpBwCluster, arUserconfigBwCluster

global ar
if(isempty(ar))
    error('please initialize by arInit')
end
if ~isfield(ar.config, 'cluster')
    arUserConfigBwCluster
end
sshLoginString = [ar.config.cluster.username '@' ar.config.cluster.loginNodeUrl];

if ischar(jobId) && strcmp(jobId, 'all')
    partition = 'multi';
    system(['ssh ',sshLoginString, ' "scancel --user=' ar.config.cluster.username ' --partition=' partition '"']);
elseif isnumeric(jobId)
    toKillStr = sprintf('%i ', jobId);
    system(['ssh ',sshLoginString, ' "scancel ' toKillStr '"']);
else
    error('Invalid JOBID. Must be numeric or ''all''')
end
end
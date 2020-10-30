function varargout = arJobStatusBwCluster(cjId, silent)

% status = arJobStatusBwCluster(cjId, [silent])
%
% Get queue overview and check status of jobs on the cluster corresponding
% to cjId
%
%   [silent]    return status [false]
%
% See also: arHelpBwCluster, arUserConfigBwCluster,
% arDownloadResultsBwCluster

if ~exist('silent') || isempty(silent)
    silent = true;
end

%% CONFIG
global ar
if(isempty(ar))
    error('please initialize by arInit')
end
if ~isfield(ar.config, 'cluster')
    arUserConfigBwCluster
end
if ~isfield(ar.config, 'cluster') || ~isfield(ar.config.cluster, 'mapping') || ~isfield(ar.config.cluster.mapping, cjId) || ~isfield(ar.config.cluster.mapping.(cjId), 'jobIds')
    error('No cjId -> jobId mapping found in ar struct')
end

sshLoginString = [ar.config.cluster.username '@' ar.config.cluster.loginNodeUrl];

%% SSH WORK
% retrieve mapping of cjId to corresponding cluster jobIds
try
    myjobs = ar.config.cluster.mapping.(cjId).jobIds;
catch
    error('Did not find mapping for provided cjId')
end

% check presence of slurm-out files via ssh
jobIdString = cellfun(@(x) sprintf('%i',x), num2cell(myjobs'), 'UniformOutput', false);
slurmOutFiles = cellfun(@(x) sprintf('slurm-%i.out',x), num2cell(myjobs'), ...
    'UniformOutput', false);
checkIfExist = strcat('if test -f \"', slurmOutFiles, '\"; then echo \"', jobIdString, ' completed\"; else echo \"', jobIdString, ' inprogress\"; fi');
checkIfExist = strjoin(checkIfExist, '; ');
%fprintf('%s@%s''s password:\n', ar.config.cluster.username, ar.config.cluster.loginNodeUrl)
[status, cmdout] = system(['ssh ',sshLoginString, ' "squeue && cd ', ar.config.cluster.wd,'; ', checkIfExist, '"'], '-echo');

%% INTERPRET RESULTS
resPos = length(regexp(cmdout, 'completed'));
resNeg = length(regexp(cmdout, 'inprogress'));
if resPos + resNeg ~= length(myjobs)
    error('Bad file count')
end

if resPos < length(myjobs)
    fprintf('\narJobStatusBwCluster: In progress. %i of %i jobs finished\n', resPos, length(myjobs))
    status = 0;
else
    fprintf('\narJobStatusBwCluster: All jobs are finished\n')
    status = 1;
end

if ~silent
    varargout{1} = status;
end
end
function arDownloadResultsBwCluster(cjId, resFolderName)

% arDownloadResultsBwCluster(cjId, [resFolderName])
%
% Start results collection on and download results of job corresponding to 
% cjId from the cluster.
%
%   resFolderName      Name of the local results folder
%
% See also arHelpBwCluster, arUserConfigBwCluster, arKillJobBwCluster

global ar

if(isempty(ar))
    error('please initialize by arInit')
end
if ~isfield(ar.config, 'cluster')
    if exist(sprintf('clConfBackup_%s.mat', cjId), 'file')
        fprintf('Did not find field ar.config.cluster. Cluster config is loaded from backup file %s', sprintf('clConfBackup_%s.mat', cjId))
        load(sprintf('clConfBackup_%s.mat', cjId), 'clConf');
        ar.config.cluster.username = clConf.username;
        ar.config.cluster.loginNodeUrl = clConf.loginNodeUrl;
        ar.config.cluster.matlabVersion = clConf.matlabVersion;
        ar.config.cluster.d2dpath = clConf.d2dpath;
        ar.config.cluster.wd = clConf.wd;
        ar.config.cluster.mapping.(cjId).jobIds = clConf.mapping.(cjId).jobIds;
    else
        error('Did not find field ar.config.cluster or backup file. Please write config by running arUserConfigBwCluster.')
    end
end

arCheckFields(ar.config.cluster, {'username', 'loginNodeUrl', 'wd', sprintf('mapping.%s.jobIds', cjId), sprintf('mapping.%s.rescolscriptfilename', cjId)}, 2);

%% CONFIG
sshLoginString = [ar.config.cluster.username '@' ar.config.cluster.loginNodeUrl];

clusterWd = ar.config.cluster.wd;
rescolscriptfilename = ar.config.cluster.mapping.(cjId).rescolscriptfilename;
name = ar.config.cluster.mapping.(cjId).name;
arMatFilePath = [clusterWd, '/', name, '/Results/', name, '/workspace.mat'];

if ~exist('resFolderName') || isempty(resFolderName)
    resFolderName = name;
end

destDir = ['Results', filesep, resFolderName];
if ~exist('Results', 'dir'); mkdir('Results'); end
if ~exist(destDir, 'dir'); mkdir(destDir); end

%% CHECK JOB STATUS

fprintf('Checking job status...\n')
status = arJobStatusBwCluster(cjId, false);
if status == 0
    error('Jobs have not yet finished.')
else
    fprintf('Proceeding results collection on cluster...\n')
end

%% COLLECTING BASH COMMANDS
% open ssh master connection
sshstring = ['SSHSOCKET=~/.ssh/' sshLoginString ' && ',...
    'ssh -M -f -N -o ControlPath=$SSHSOCKET ' sshLoginString];

% trigger results collection
sshstring = [sshstring, ' && ',...
    'ssh -o ControlPath=$SSHSOCKET ' sshLoginString ' "cd ' clusterWd '/' name ' && module load math/matlab/' ar.config.cluster.matlabVersion ' && pwd && matlab -r \"' rescolscriptfilename(1:(end-2)) ';exit' '\""'];

% download results
sshstring = [sshstring, ' && ',...
    'scp -r -o ControlPath=$SSHSOCKET ', sshLoginString, ':', arMatFilePath, ' ', destDir, '/workspace.mat'];

% close ssh master connection
sshstring = [sshstring, ' && ', ...
    'ssh -S $SSHSOCKET -O exit ', sshLoginString];

%% EXECUTE BASH COMMANDS
[status, cmdout] = system(sshstring, '-echo');
if status == 0
    fprintf('arDownloadResultsBwCluster: Download successfull. Type ''arLoad(''%s'')'' to load the workspace. \nar.p is not automaticcaly updated with results from multistart optimization!\n', resFolderName)
else
    fprintf(2, 'arDownloadResultsBwCluster: Download failed! Check bash output above!\n')
end
end
function jobIds = arSendJobBwCluster(cjId, functionString, name, uploadFiles)

% jobIds = arSendJobBwCluster(cjId, functionString, [name], [uploadFiles])
%
% Send computing job (arFitLhsBwCluster or pleBwCluster) to cluster.
%   cjId                identifier for computing job, (e.g. 'run01')
%   name                folder name on cluster within ar.config.cluster.wd
%                       [name of the current dir]
%   functionString      function to execute on login node as string 
%                       (e.g. 'arFitLhsBwCluster(30,5)'
%   uploadFiles         upload files to cluster (arUploadToBwCluster) [false]
%
% See also: arHelpBwCluster, arUserConfigBwCluster, arUploadBwCluster, 
% arFitLhsBwCluster, pleBwCluster, arJobStatusBwCluster,
% arDownloadResultsBwCluster, arKillJobBwCluster

global ar
if ~exist('name') || isempty(name)
    [~,name,~]=fileparts(pwd);
end

if ~exist('cjId') || isempty(cjId)
    cjId = datestr(now,30);
end
if(isempty(ar))
    error('please initialize by arInit')
end
if ~isfield(ar.config, 'cluster')
    arUserConfigBwCluster
end

if isempty(regexp(functionString, '(arFitLhsBwCluster|pleBwCluster)\(.*\)', 'ONCE'))
    error('Invalid function string. Must be of the form ''arFitLhsBwCluster(...)'' or ''pleBwCluster(...)''')
end
if ~exist('uploadFiles') || isempty(uploadFiles)
    uploadFiles = false;
end

if uploadFiles
    arUploadBwCluster(name)
end

%% CONFIG
sshLoginString = [ar.config.cluster.username '@' ar.config.cluster.loginNodeUrl];
scriptfilename = ['m_', cjId, '_', name, '_clusterInit.m'];
rescolscriptfilename = ['m_', cjId, '_', name, '_resultsCollection.m'];
clusterWd = ar.config.cluster.wd;

ar.config.cluster.mapping.(cjId).rescolscriptfilename = rescolscriptfilename;
ar.config.cluster.mapping.(cjId).name = name;


%% COLLECTING BASH COMMANDS
% open ssh master connection
sshstring = ['SSHSOCKET=~/.ssh/' sshLoginString ' && ',...
    'ssh -M -f -N -o ControlPath=$SSHSOCKET ' sshLoginString];

% generate and upload matlab script
writeClusterInitFile(ar.config.cluster.d2dpath, name, scriptfilename, rescolscriptfilename, functionString);
sshstring = [sshstring, ' && ',...
    'scp -o ControlPath=$SSHSOCKET ', scriptfilename, ' ', sshLoginString, ':', clusterWd, '/', name];

% run matlab script
sshstring = [sshstring, ' && ',...
    'ssh -o ControlPath=$SSHSOCKET ' sshLoginString ' "cd ' clusterWd '/' name ' && module load math/matlab/' ar.config.cluster.matlabVersion ' && pwd && matlab -r \"' scriptfilename(1:(end-2)) '\""'];

% close ssh master connection
sshstring = [sshstring, ' && ', ...
    'ssh -S $SSHSOCKET -O exit ', sshLoginString];

%% Execute bash commands
[status, cmdout] = system(sshstring, '-echo');
if status == 0
    fprintf('arSendJobToBwCluster: All jobs submitted\n')
else 
    fprintf(2, 'arSendJobToBwCluster: Job submission failed. Check bash output above!\n')
end

%% Get jobIds from output
jobIds = str2double(regexp(cmdout, '(?<=Submitted batch job )[0-9]*', 'match'));
ar.config.cluster.mapping.(cjId).jobIds = jobIds;
ar.config.cluster.mapping.(cjId).name = name;

clConf = ar.config.cluster;
save(sprintf('clConfBackup_%s.mat', cjId), 'clConf')
end

function writeClusterInitFile(d2dpath, name, scriptfilename, rescolscriptfilename, functionString)
mcode = {
 %   ['cd ',clusterWd, '/', name], ...
    ['addpath(''',d2dpath,''');'], ...
    'arInit;', ...
    ['arLoad(''',name,''');'],...
    '', ...
    'if ~isfile([ar.fkt ''.mexa64'']); arRecompile; end',...
    ['resColFun = ' functionString ';'], ...
    '',...
    'x = strsplit(resColFun, ''/'');',...
    'x = strsplit(x{end},''.'');',...
    ['save(''' rescolscriptfilename ''',''resColFun'');'],...
    ['fileID = fopen(''' rescolscriptfilename ''', ''w'');'],...
    ['fprintf(fileID, x{1});'],...
    ['fclose(fileID);'],...
    'exit',...
    };

%filename = confCl.file_init;
fid = fopen(scriptfilename,'w');
for i=1:length(mcode)
    fprintf(fid,'%s\n',mcode{i});
end
end

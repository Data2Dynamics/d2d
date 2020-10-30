function arUploadBwCluster(name)

% arUploadBwCluster([name])
%
% Upload current workspace with all necessary files for recompilation to
% cluster. 
%
%   name        Name of the results folder
%
% If executed on Linux systems, arSimuCalcFun will also be uploaded to the cluster to 
% remove the need for recompilation.
%
% See also arHelpBwCluster, arUserConfigBwCluster, arSendJobBwCluster

global ar
if(isempty(ar))
    error('please initialize by arInit')
end
if ~isfield(ar.config, 'cluster')
    arUserConfigBwCluster
end
if ~exist('name') || isempty(name)
    name2 = strsplit(ar.config.savepath, filesep);
    name = [name2{end}, datestr(now,30)];
end


%% CONFIG
confCl.name = name;
confCl.arMatFileName = 'workspace.mat';

confCl.arPathLoc = [confCl.arMatFileName];      % path to ar struct on local

confCl.wd = ar.config.cluster.wd; % pwd on cluster
confCl.arPathCl = [confCl.wd, '/', confCl.arMatFileName];        % path to ar struct on cluster

sshLoginString = [ar.config.cluster.username '@' ar.config.cluster.loginNodeUrl];


%% COLLECTING BASH COMMANDS
% open ssh master connection
sshstring = ['SSHSOCKET=~/.ssh/' sshLoginString ' && ',...
    'ssh -M -f -N -o ControlPath=$SSHSOCKET ' sshLoginString];

% model and data files from backup dirs
sshstring = [sshstring, ' && ', ...
    getSSHstrings(ar.setup.modelfiles, ar.setup.backup_model_folder, ...
    ar.setup.backup_model_folder_local, ...
    [confCl.wd, '/Models/Backup/'], sshLoginString)];
sshstring = [sshstring, ' && ',...
    getSSHstrings(ar.setup.datafiles, ar.setup.backup_data_folder, ...
    ar.setup.backup_data_folder_local, ...
    [confCl.wd, '/Data/Backup/'], sshLoginString)];

% copy simucalc fun if on linux
if isunix && ~ismac
sshstring = [sshstring, ' && ',...
    'scp -o ControlPath=$SSHSOCKET ', ar.fkt, '.mexa64 ', sshLoginString, ':', confCl.wd];
end

% ar struct
save(confCl.arPathLoc, 'ar')
sshstring = [sshstring, ' && ',...
    'ssh -o ControlPath=$SSHSOCKET ', sshLoginString, ' "cd ', confCl.wd, ' && if ! test -d \"Results\"; then mkdir Results; fi" && ', ...
    'ssh -o ControlPath=$SSHSOCKET ', sshLoginString, ' "cd ', confCl.wd, ' && mkdir Results/', confCl.name, '" && ', ...
    'scp -o ControlPath=$SSHSOCKET ', confCl.arPathLoc, ' ', sshLoginString, ':', confCl.wd, '/Results/', confCl.name];

% close ssh master connection
sshstring = [sshstring, ' && ', ...
    'ssh -S $SSHSOCKET -O exit ', sshLoginString];


%% EXECUTE BASH COMMANDS
[status, cmdout] = system(sshstring, '-echo');
if status == 0
    fprintf('arUploadToBwCluster: Upload successfull\n')
else
    fprintf(2, 'arUploadToBwCluster: Upload failed/incomplete! Check bash output above!\n')
end

end

function sshstring = getSSHstrings(setupFiles, setupBackupFolder, setupBackupFolderLocal, destination, sshLoginString)
% copy files from ar.setup to destination
%   setupFiles: ar.setup.modelfiles
%   setupBackupFolder: ar.setup.backup_model_folder
%   setupBackupFolderLocal: ar.setup_backup_folder_local
sshstring = 'TRUE'; %dummy

for io = 1:length(setupFiles)
    if ~isempty(setupFiles{io})
        if iscell(setupFiles{io})
            for ii = 1:length(setupFiles{io})
                from = setupBackupFolderLocal{io}{ii};
                to = [destination];
                
                sshstring = [sshstring, ' && ',...
                    'ssh -o ControlPath=$SSHSOCKET ' sshLoginString ' "mkdir -p ', to, '"',...
                    ' && ',...
                    'scp -r -o ControlPath=$SSHSOCKET ', from, ' ', sshLoginString, ':', to];
            end
        else
            from = setupBackupFolder{io};
            to = [destination];
            
            sshstring = [sshstring, ' && ',...
                'ssh -o ControlPath=$SSHSOCKET ', sshLoginString, ' "mkdir -p ', to, '"',...
                ' && ',...
                'scp -r -o ControlPath=$SSHSOCKET ', from, ' ', sshLoginString, ':', to];
        end
    end
end
end
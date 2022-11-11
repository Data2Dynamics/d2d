% arSave(name, [withSyms], [useDateString])
%
% arSave saves the current workspace, i.e. model struct and last ple results
% to repository, or return base path of last arSave.
%
%   arSave
%   arSave(name)
%   arSave(name, withSyms)
%   basepath = arSave(--)
%
%   name            char indicating the directory name where to save the
%                   workspace - if left empty or set to 'current', current path
%                   in ar.config.savepath is used
%   withSyms        optional logical indicating wheter symbols in the ar struct
%                   shall be saved or not. withSyms true needs the symbolic
%                   toolbox. Inclusion of symbols can be memory intensive [false]
%   useDateString   prepend unique datestring to results directory name [true]
%
%   basepath  (optional) base path of last arSave to name
%
% Example
%    Load ABC model and data, compile model and save ar struct to
%    repository:
%    arLoadModel('ABC_model');
%    arLoadData('ABC_data_BCobs');
%    arCompileAll();
%    Save workspace to the current repository with path in ar.config.savepath%
%    arSave('current')
%    Save workspace to the new repository with name 'ABC_test'
%    arSave('ABC_test')
%
% See also arSaveParOnly

function basepath = arSave(name, withSyms, useDateString)

global ar

% Invalidate the cache going into the save
arCheckCache(1);

% update checkstrs:
ar = arUpdateCheckstr(ar, true);

if(~exist('withSyms','var') || isempty(withSyms))
    withSyms = false;
end
if ~exist('useDateString') || isempty(useDateString)
    useDateString = true;
end


if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end



if(isempty(ar.config.savepath)) % never saved before, ask for name
    if(~exist('name','var'))
        name = input('enter new repository name addition: ', 's');
    end
    if(~isempty(name))
        if useDateString
            ar.config.savepath = ['./Results/' datestr(now, 30) '_' name];
        else
            ar.config.savepath = ['./Results/' name];
        end
    else
        ar.config.savepath = ['./Results/' datestr(now, 30) '_noname'];
    end
    
    if isfield(ar.config,'backup_modelAndData') && ar.config.backup_modelAndData
        createBackupFiles
    end
    
    if (~ar.config.lightSave)
        arSaveFull(ar,withSyms);
    else
        warning( 'ar.config.lightSave is set to true. Only saving parameter set' );
    end
    arSaveParOnly(ar, ar.config.savepath);
else
    if(exist('name','var')) % saved before, but new name give
        if(~strcmpi(strtrim(name), 'current'))
            if(~isempty(name))
                if useDateString
                    ar.config.savepath = ['./Results/' datestr(now, 30) '_' name];
                else
                    ar.config.savepath = ['./Results/' name];
                end
            else
                ar.config.savepath = ['./Results/' datestr(now, 30) '_noname'];
            end
        end
        
        if isfield(ar.config,'backup_modelAndData') && ar.config.backup_modelAndData
            createBackupFiles
        end
        
        if (~ar.config.lightSave)
            arSaveFull(ar,withSyms);
        else
            warning( 'ar.config.lightSave is set to true. Only saving parameter set' );
        end
        arSaveParOnly(ar, ar.config.savepath);
        
    else
        if(nargout == 0) % saved before, ask for new name give
            name = input(sprintf('enter new repository name addition [%s]: ', ...
                ar.config.savepath), 's');
            
            if(~isempty(name))
                ar.config.savepath = ['./Results/' datestr(now, 30) '_' name];
            end
            
            if isfield(ar.config,'backup_modelAndData') && ar.config.backup_modelAndData
                createBackupFiles
            end
            
            if (~ar.config.lightSave)
                arSaveFull(ar,withSyms);
            else
                warning( 'ar.config.lightSave is set to true. Only saving parameter set' );
            end
            arSaveParOnly(ar, ar.config.savepath);
        else
            if(~exist(ar.config.savepath, 'dir'))
                % saved before, path output requested,
                % however save path does not exist anymore
                
                if isfield(ar.config,'backup_modelAndData') && ar.config.backup_modelAndData
                    createBackupFiles
                end
                
                if (~ar.config.lightSave)
                    arSaveFull(ar,withSyms);
                else
                    warning( 'ar.config.lightSave is set to true. Only saving parameter set' );
                end
                arSaveParOnly(ar, ar.config.savepath);
            end
        end
    end
end

if(nargout>0)
    basepath = ar.config.savepath;
end




function createBackupFiles
% Create backup files (data and model) for arRecompile
%   - For new first saves after compiling copy files from Data/ and Models/
%   - For resaved workspace copy backup files from previous workspace

global ar

% Create folders if not yet done
if(~exist(ar.config.savepath, 'dir'))
    mkdir(ar.config.savepath)
end

if ~isfolder([ar.config.savepath '/Backup'])
    mkdir([ar.config.savepath '/Backup']);
end
if ~isfolder([ar.config.savepath '/Backup/' ar.checkstr])
    mkdir([ar.config.savepath '/Backup/' ar.checkstr]);
end

dataBackupDir = [ar.config.savepath '/Backup/' ar.checkstr '/Data'];
modelBackupDir = [ar.config.savepath '/Backup/' ar.checkstr '/Models'];
if ~isfolder(dataBackupDir)
    mkdir(dataBackupDir);
end
if ~isfolder(modelBackupDir)
    mkdir(modelBackupDir);
end




% Check if this is first save after compilation, 
if isfield(ar,'setup') % not available for old workspaces
    if ~isfield(ar.setup,'backup_data_folder')
        % copy files from Data/ and Models/
        arFprintf(2,'Creating backup.')
        
        ar.setup.backup_data_folder = cell(size(ar.setup.datafiles));
        ar.setup.backup_model_folder = cell(size(ar.setup.modelfiles));
        
        for i=1:length(ar.setup.datafiles)
            for j=1:length(ar.setup.datafiles{i})
                %             old_backup_data_folder = ar.setup.backup_data_folder{i}{j};
                ar.setup.backup_data_folder{i}{j} = fullfile(pwd,[dataBackupDir, '/']);% fullfile to prevent mixing of \ and /
                ar.setup.backup_data_folder_local{i}{j} = [dataBackupDir, '/'];
                [~,file,ext] = fileparts(ar.setup.datafiles{i}{j});
                source = ar.setup.datafiles{i}{j};
                %             source = [old_backup_data_folder,file,ext];
                target = [ar.setup.backup_data_folder{i}{j},file,ext];
                if ~isempty(source) && ~isempty(target) && strcmp(fullfile(strrep(source,pwd,'.')),fullfile(strrep(target,pwd,'.')))~=1
                    try
                        copyfile(source,target);
                    catch
                        arFprintf(2,'Error while copying files \n from %s\n to %s\n',source,target)
                        error('Failed to create backup files. Turning off ar.config.backup_modelAndData might be an option.')
                    end
                end
                
            end
        end
        
        for i=1:length(ar.setup.modelfiles)
            if ~isempty(ar.setup.modelfiles{i})
                %             old_backup_model_folder = ar.setup.backup_model_folder{i};
                ar.setup.backup_model_folder{i} = fullfile(pwd,[modelBackupDir, '/']);% fullfile to prevent mixing of \ and /
                ar.setup.backup_model_folder_local{i} = [modelBackupDir, '/'];
                [~,file,ext] = fileparts(ar.setup.modelfiles{i});
                source = ar.setup.modelfiles{i};
                %             source = [old_backup_model_folder,file,ext];
                target = [ar.setup.backup_model_folder{i},file,ext];
                if ~isempty(source) && ~isempty(target) && strcmp(fullfile(strrep(source,pwd,'.')),fullfile(strrep(target,pwd,'.')))~=1
                    try
                        copyfile(source,target);
                    catch
                        arFprintf(2,'Error while copying files \n from %s\n to %s\n',source,target)
                        error('Failed to create backup files. Turning off ar.config.backup_modelAndData might be an option.')
                    end
                end
            end
        end
        
        % If this is not first save, and ar.setup.backup_data_folder is not empty
    elseif ~isempty(ar.setup.backup_data_folder)
        % copy files from previous workspace
        old_backup_data_folder = ar.setup.backup_data_folder_local;
        old_backup_model_folder = ar.setup.backup_model_folder_local;
        
        
        for i=1:length(ar.setup.datafiles)
            for j=1:length(ar.setup.datafiles{i})
                %             old_backup_data_folder = ar.setup.backup_data_folder{i}{j};
                ar.setup.backup_data_folder{i}{j} = fullfile(pwd,[dataBackupDir, '/']);% fullfile to prevent mixing of \ and /
                ar.setup.backup_data_folder_local{i}{j} = [dataBackupDir, '/'];
                [~,file,ext] = fileparts(ar.setup.datafiles{i}{j});
                source = [old_backup_data_folder{i}{j},file,ext];
                target = [ar.setup.backup_data_folder{i}{j},file,ext];
                if ~isempty(source) && ~isempty(target) && strcmp(fullfile(strrep(source,pwd,'.')),fullfile(strrep(target,pwd,'.')))~=1
                    try
                        [suc,mess,messID]=copyfile(source,target);
                    catch
                        source
                        target
                        mess
                        arFprintf(2,'Error while copying files \n from %s\n to %s\n',source,target)
                        error('Failed to create backup files. Turning off ar.config.backup_modelAndData might be an option.')
                    end
                end
            end
        end
        
        for i=1:length(ar.setup.modelfiles)
            if ~isempty(ar.setup.modelfiles{i})
                %             old_backup_model_folder = ar.setup.backup_model_folder{i};
                ar.setup.backup_model_folder{i} = fullfile(pwd,[modelBackupDir, '/']);% fullfile to prevent mixing of \ and /
                ar.setup.backup_model_folder_local{i} = [modelBackupDir, '/'];
                [~,file,ext] = fileparts(ar.setup.modelfiles{i});
                source = [old_backup_model_folder{i},file,ext];
                target = [ar.setup.backup_model_folder{i},file,ext];
                if ~isempty(source) && ~isempty(target) && strcmp(fullfile(strrep(source,pwd,'.')),fullfile(strrep(target,pwd,'.')))~=1
                    try
                        copyfile(source,target);
                    catch
                        arFprintf(2,'Error while copying files \n from %s\n to %s\n',source,target)
                        error('Failed to create backup files. Turning off ar.config.backup_modelAndData might be an option.')
                    end
                end
            end
        end
        
    else
        error('ar.setup.backup_data_folder empty. Creating backup files was not possible. Try turning off ar.config.backup_modelAndData.')
    end
    
end




% full save
function arSaveFull(ar,withSyms)
% global ar

% remove storage-consuming fields in global ar
% that are non-essential
arCompress;

% check if dir exists
if(~exist(ar.config.savepath, 'dir'))
    mkdir(ar.config.savepath)
end

% remove symy for saving
if(~withSyms)
    for jm = 1:length(ar.model)
        ar.model(jm).sym = [];
        for jc = 1:length(ar.model(jm).condition)
            ar.model(jm).condition(jc).sym = [];
        end
        if(isfield(ar.model(jm), 'data'))
            for jd = 1:length(ar.model(jm).data)
                ar.model(jm).data(jd).sym = [];
            end
        end
    end
end
arDoSaving(ar)

% Copy mex file to the save folder
if ~isfield(ar.config, 'saveMexFile')
    ar.config.saveMexFile = true;
end
if ar.config.saveMexFile
    % Read out file ending of mex file for current os
    if ismac
        osext = 'mexmaci64';
    elseif isunix
        osext = 'mexa64';
    elseif ispc
        osext = 'mexw64';
    else
        disp('Platform not supported')
    end
    fkt = which([ar.fkt '.' osext]);
    [~, fkt_name, fkt_extension] = fileparts(fkt);
    
    if ~isempty(fkt) && ~exist([ar.config.savepath '/' fkt_name fkt_extension],'file')
        copyfile(fkt , ar.config.savepath)
    end
end



% This function does not use ar as global variables for
% being able to do manipulations without altering the global variables.
function arDoSaving(ar)
ar = arDeleteGraphicsHandles(ar);

% warning off MATLAB:Figure:FigureSavedToMATFile
w = whos('ar');
if w.bytes>1e9;
    save([ar.config.savepath '/workspace.mat'],'ar','-v7.3'); % can be 100x larger than -v7 (e.g. for the Chen model)
else
    save([ar.config.savepath '/workspace.mat'],'ar','-v7');
end
% warning on MATLAB:Figure:FigureSavedToMATFile

fprintf('workspace saved to file %s\n', [ar.config.savepath '/workspace.mat']);








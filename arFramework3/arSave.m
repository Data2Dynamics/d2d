% arSave(name, [withSyms])
% 
% arSave saves the current workspace, i.e. model struct and last ple results
% to repository, or return base path of last arSave.
%
%   arSave
%   arSave(name)
%   arSave(name, withSyms)
%   basepath = arSave(--)
%
%   name      char indicating the directory name where to save the
%             workspace - if left empty or set to 'current', current path 
%             in ar.config.savepath is used
%   withSyms  optional logical indicating wheter symbols in the ar struct 
%             shall be saved or not. withSyms true needs the symbolic 
%             toolbox. Inclusion of symbols can be memory intensive [false]
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

function basepath = arSave(name, withSyms)

global ar

% Invalidate the cache going into the save
arCheckCache(1);

% update checkstrs:
ar = arUpdateCheckstr(ar, true);

if(~exist('withSyms','var'))
    withSyms = false;
end
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end


if(isempty(ar.config.savepath)) % never saved before, ask for name
    if(~exist('name','var'))
        name = input('enter new repository name addition: ', 's');
    end
    if(~isempty(name))
        ar.config.savepath = ['./Results/' datestr(now, 30) '_' name];
    else
        ar.config.savepath = ['./Results/' datestr(now, 30) '_noname'];
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
                ar.config.savepath = ['./Results/' datestr(now, 30) '_' name];
            else
                ar.config.savepath = ['./Results/' datestr(now, 30) '_noname'];
            end
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

% full save
function arSaveFull(ar,withSyms)
% global ar

% remove storage-consuming fields in global ar
% that are non-essential
arCompress;

% check is dir exists
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
    fkt = which(ar.fkt);
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








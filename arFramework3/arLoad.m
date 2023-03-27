% arLoad
% arLoad(workspace_name)
% arLoad(workspace_name,doAutomaticRecompile, sortModus)
% 
% arLoad loads a model struct with latest ple results from the repository
% 
%   workspace_name    name of the result folder or
%                     ar.config.savepath or
%                     an integer selected from the list displayed by 
%                     call arLoad
% 
%   doAutomaticRecompile [0] 
%         If ar.fkt is not found (e.g. due to compiling
%         previously on another OS, arRecompile can be started
%         automatically if set to 1
% 
% sortModus      sorting of the workspaces (passed to fileChooser.m > fileList.m)
%               'none' [default]
%               'chi2'
%               'checkstr'
%
% Examples
%    Load workspace from current/last path at ar.config.savepath
%    arLoad
%    arLoad(1)
%    arLoad(ar.config.savepath)
%    Load a struct from an existing folder
%    arLoad('NameOfAnExistingResultFolder')
%
% arLoad([],[],'chi2')
% arLoad([],[],'checkstr')
%
% See also arLoadLatest, arLoadData, arLoadFilename 

function arLoad(workspace_name,doAutomaticRecompile,sortModus)
if ~exist('sortModus','var') || isempty(sortModus)
    sortModus = 'none';
end
if ~exist('doAutomaticRecompile','var') || isempty(doAutomaticRecompile)
    doAutomaticRecompile = 0;
end

if ~exist('Results','dir')
    error('No results folder exist. arLoad can only be executed in a D2D working directory.')
end
if exist('arCheck','file')==0
    arInit;
end
arCheck;

global ar


% set the two variables also as global in the command line workspace:
evalin('base','clear ar');  
evalin('base','global ar');  

if(~exist('workspace_name', 'var') || isempty(workspace_name))
    [~, workspace_name] = fileChooser('./Results', 1, true, [], sortModus);
elseif(isnumeric(workspace_name)) % workspace_name is the file-number
    [~, ~, file_list] = fileChooser('./Results', 1, -1, [], sortModus);    
    workspace_name = file_list{workspace_name};
elseif(ischar(workspace_name)) 
    [~,workspace_name]=fileparts(workspace_name);    % remove path
end
if contains(workspace_name,'*')
    [~, workspace_name] = fileChooser('Results', 1, true, replace(workspace_name,'*',''), sortModus);    
end

Stmpload = load(['./Results/' workspace_name '/workspace.mat']);
ar = Stmpload.ar;

if(strcmp(ar.config.savepath,['./Results/' workspace_name])~=1)
    ar.config.savepath = ['./Results/' workspace_name];
end

if ~isfield(ar,'ple') || isempty(ar.ple) || ~isfield(ar.ple,'ps')
    ar.ple = pleLoad(ar); % if an old/deprecated workspace is available
end

fprintf('workspace loaded from file %s\n', workspace_name);

% Check if the mex file is in the working directory, if not try to copy it
% from the savefolder
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

if isempty(fkt)
    fprintf([ar.fkt '.' osext ' not found. Try to copy it from the savefolder. \n'])  
    files = dir([ar.config.savepath '/' ar.fkt '.' osext]);
    cf_succeed = 0;
    for idf = 1:length(files)
        copyfile([ar.config.savepath '/' files.name ] , pwd)
        fprintf(' done.\n')
        cf_succeed = 1;
    end
    if ~cf_succeed
        if doAutomaticRecompile
            arRecompile
        else
            fprintf([ar.fkt '.' osext '  NOT found. Please compile via Setup or ''arRecompile'' or set option doAutomaticRecompile=1.\n'])
        end
    end
end

% Make sure we have all the necessary fields
ar = arInitFields(ar);




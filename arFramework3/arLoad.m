% arLoad(workspace_name)
% 
% load model struct and last ple results
% 
%   workspace_name    either 
%                 - the name of the result folder 
%                 - ar.config.savepath
%                 - or the number as displayed by arLoad without argument
%   Examples:
% arLoad
% arLoad(1)
% arLoad(ar.config.savepath)
% arLoad('NameOfAnExistingResultFolder')

function arLoad(workspace_name)
if exist('arCheck','file')==0
    arInit;
end
arCheck;

global ar

% set the two variables also as global in the command line workspace:
evalin('base','clear ar');  
evalin('base','global ar');  

if(~exist('workspace_name', 'var') || isempty(workspace_name))
    [~, workspace_name] = fileChooser('./Results', 1, true);
elseif(isnumeric(workspace_name)) % workspace_name is the file-number
    [~, ~, file_list] = fileChooser('./Results', 1, -1);    
    workspace_name = file_list{workspace_name};
elseif(ischar(workspace_name)) 
    [~,workspace_name]=fileparts(workspace_name);    % remove path
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
fkt = which(ar.fkt);

if isempty(fkt)
    disp([ar.fkt ' not found. Try to copy it from the savefolder...'])
    files = dir([ar.config.savepath '/' ar.fkt '.*']);
    cf_succeed = 0;
    for idf = 1:length(files)
        copyfile([ar.config.savepath '/' files.name ] , pwd)
        disp(['Copied ' files.name ' from ' ar.config.savepath ' to the working directory.'])
        cf_succeed = 1;
    end
    if ~cf_succeed
        disp([ar.fkt ' not found in ' ar.config.savepath])
    end
end

% Make sure we have all the necessary fields
ar=arInitFields(ar);




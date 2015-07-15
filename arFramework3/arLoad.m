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
arCheck;

global ar
global pleGlobals

% set the two variables also as global in the command line workspace:
evalin('base','clear ar pleGlobals');  
evalin('base','global ar pleGlobals');  

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
    ar.config.savepath = ['./Results/' workspace_name]
end
    

fprintf('workspace loaded from file %s\n', workspace_name);

try
    pleGlobals = pleLoad(ar);
catch
    fprintf(1,'No valid PLE workspace found!\n');
    clear pleGlobals;
end





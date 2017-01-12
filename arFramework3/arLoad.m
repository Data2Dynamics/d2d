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
if ~isfield(ar.config,'add_c')  % if ar has been created with an old D2D version
    ar.config.add_c = 50;
end

% new:
if(isfield(Stmpload,'pleGlobals'))
    pleGlobals = Stmpload.pleGlobals;  % is overwritten, if ple in PLE/result.mat is available and finished.
end
% end new.
if(strcmp(ar.config.savepath,['./Results/' workspace_name])~=1)
    ar.config.savepath = ['./Results/' workspace_name];
end

fprintf('workspace loaded from file %s\n', workspace_name);

try
    ple = pleLoad(ar);
    % as before:
    if isempty(pleGlobals)  %no pleGlobals in ar.config.savepath/workspace.mat, use the PLEs in the ar.config.savepath/results.mat
        pleGlobals = ple;
    elseif(~isempty(ple) && ple.finished) % A finished calculation in ar.config.savepath/results.mat is available, it has priority overwrites potential ples in ar.config.savepath/workspace.mat
        pleGlobals = ple;
    % new:
    else % pleGlobals in ar.config.savepath/results.mat available, but not finished: Keep pleGlobals as saved by arSave, i.e. use/keep ples from ar.config.savepath/workspace.mat
        fprintf('Calculation of temporary PLE in %s/PLE/result.mat is not finished.\n>Load PLE from %s/workspace.mat instead.\n',workspace_name,workspace_name);
    end
    % end new
catch
    fprintf(1,'No valid PLE workspace found!\n');
    clear pleGlobals;
end

% Make sure we have all the necessary fields
ar=arInitFields(ar);




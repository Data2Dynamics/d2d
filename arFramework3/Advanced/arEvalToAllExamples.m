% arEvalToAllExamples(fun)
% 
% arEvalToAllExamples(fun, [], [], [], varargin)
% 
% arEvalToAllExamples(fun, result_suffix, workspace_pattern, depth, varargin)
% 
% 
%   Loads recursively the latest result workspaces and applies the function
%   fun. 
% 
%   The function checks whether a folder is a D2D working folder by
%   checking existance of the folders 'Results' and 'Models'.
% 
% 
% 
%   fun             a D2D or user-defined function 
% 
%   result_suffix   used as argument for arSave after fun is finished
% 
%               If result_suffix = 'current' then the loaded workspace is
%               used
% 
%   workspace_pattern   used as argument for arLoadLatest to load a
%                       workspace
% 
%   depth       folder depth for which the fucntion should work.
%               The argument depth is given to list_files_recursive(depth).
% 
%               depth=1 applies fun in the current folder and in all
%               existing subfolders (but not in subsubfolder)
% 
%               depth = Inf means working recursively through all folders
%               (default)
% 
%   vargin      additional arguments passed to fun
% 
% 
% Example:
% arEvalToAllExamples('arFitLHS','LHS100',[],1,100);


function arEvalToAllExamples(fun, result_suffix, workspace_pattern, depth, varargin)
if ~exist('depth','var') || isempty(depth)
    depth = Inf;
end
if ~exist('result_suffix','var') || isempty(result_suffix)
    result_suffix = 'arEvalToAllExamples';
end
if ~exist('workspace_pattern','var') || isempty(workspace_pattern)
    workspace_pattern = '';
end

[~,folders]=list_files_recursive(depth);

pw = pwd;

for i=1:length(folders)
    try
        cd(folders{i});
        fprintf('Folder %s ... \n',folders{i});
        
        if isdir('Results') && isdir('Models') % if this seems to be a D2D modelling folder;
            applyCore(fun,result_suffix,workspace_pattern,varargin{:});
        end
        
    catch ERR
        cd(pw)
        rethrow(ERR)
    end
end


function applyCore(fun,result_suffix,workspace_pattern,varargin)
global ar
status = arLoadLatest(workspace_pattern);
if status
    fprintf('Evaluating %s for workspace %s ... \n',char(fun),ar.config.savepath);
    feval(fun,varargin{:})
    arSave(result_suffix)
end

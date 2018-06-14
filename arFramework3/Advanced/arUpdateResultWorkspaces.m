%   arUpdateResultWorkspaces
% 
% This function loads all available workspaces in ./Results and saves it to
% update information saved to the latest version. The old workspaces are
% overwritten.
% 
% As an example, the variables saved in the workspace
% 'workspace_pars_only.mat' might change over time.
% 


function arUpdateResultWorkspaces

global ar
arIn = arDeepCopy(ar);

[~, ~, file_list] = fileChooser('./Results', [], -1);

try
    for i=1:length(file_list)
        fprintf('Processing %s ...\n',file_list{i});
        if exist(['Results',filesep,file_list{i},filesep,'workspace.mat'],'file')
            arLoad(file_list{i});
            if ~isempty(ar.model) && isfield(ar.model(1),'sym')
                withsyms = ~isempty(ar.model(1).sym);
            else
                withsyms = [];
            end
            arSave('current',withsyms);
        elseif exist(['Results',filesep,file_list{i},filesep,'workspace_pars_only.mat'],'file')
            tmp = load(['Results',filesep,file_list{i},filesep,'workspace_pars_only.mat']);
            arSaveParOnly(tmp.ar,['Results',filesep,file_list{i}]);
        else
            warning('No D2D workspace found in folder Results%s%s.',filesep,file_list{i});
        end
    end
catch ERR
    ar = arIn;
    rethrow(ERR)
end

ar = arIn;

    
    
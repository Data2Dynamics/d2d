%   arUpdateResultWorkspaces
% 
% This function loads all available workspaces in ./Results and saves it to
% update information saved to the latest version. The old workspaces are
% overwritten.
% 
% As an example, the variables saved in the workspace
% 'workspace_pars_only.mat' might change over time.
% 
% After the function is applied, the time stamps of the files change,
% but the order is preserved, i.e. arLoadLatest loads the same workspaces.


function arUpdateResultWorkspaces

global ar
arIn = arDeepCopy(ar);

[~, ~, file_list] = fileChooser('./Results', [], -1);

% determine the date-order of the folders in order to be sure that the
% order is preserved
dates = NaN(size(file_list));
for i=1:length(file_list)
    d = dir(['Results',filesep,file_list{i}]);
    names = {d.name};
    ind = strmatch('.',names,'exact');
    dates(i) = datenum(d(ind).date);
    [~,rf] = sort(dates); % sorting according to date  
end

try
    for i=1:length(file_list)
        fprintf('Processing %s ...\n',file_list{rf(i)});
        if exist(['Results',filesep,file_list{rf(i)},filesep,'workspace.mat'],'file')
            arLoad(file_list{rf(i)});
            if ~isempty(ar.model) && isfield(ar.model(1),'sym')
                withsyms = ~isempty(ar.model(1).sym);
            else
                withsyms = [];
            end
            arSave('current',withsyms);
        elseif exist(['Results',filesep,file_list{rf(i)},filesep,'workspace_pars_only.mat'],'file')
            tmp = load(['Results',filesep,file_list{rf(i)},filesep,'workspace_pars_only.mat']);
            arSaveParOnly(tmp.ar,['Results',filesep,file_list{rf(i)}]);
        else
            warning('No D2D workspace found in folder Results%s%s.',filesep,file_list{rf(i)});
        end
    end
catch ERR
    ar = arIn;
    rethrow(ERR)
end

ar = arIn;

    
    
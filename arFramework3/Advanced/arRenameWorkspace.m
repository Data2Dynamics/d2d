% Rename an existing saved workspace
%
% arRenameWorkspace(old_name,new_name)

function arRenameWorkspace(old_name,new_name)

global ar
global ar_tmp
ar_tmp = ar;

if(~exist('old_name','var'))
        [~, old_name] = fileChooser('./Results', 1, true);
end

Stmpload = load(['./Results/' old_name '/workspace.mat']);
ar = Stmpload.ar;

if(~exist('new_name','var'))
        new_name = input('enter new repository name addition: ', 's');
end

name = strrep(old_name,old_name(17:end),new_name);

if ~strcmpi(old_name,name)
    new_savepath = ['./Results/' name];
    move = input(sprintf('do you really want to rename the workspace \"%s\" into \"%s\"?\nY/[N] ', old_name, name),'s');
    
    if strcmpi(move,'y')
        movefile(ar.config.savepath, new_savepath);
        ar.config.savepath = new_savepath;
        save([ar.config.savepath '/workspace.mat'],'ar','-v7.3');
        
        % save only parameters
        ar2 = struct([]);
        ar2(1).pLabel = ar.pLabel;
        ar2.p = ar.p;
        ar2.qLog10 = ar.qLog10;
        ar2.qFit = ar.qFit;
        ar2.lb = ar.lb;
        ar2.ub = ar.ub;
        ar2.type = ar.type;
        ar2.mean = ar.mean;
        ar2.std = ar.std;
        arSaveParOnly(ar2, ar.config.savepath);
    end
end

ar = ar_tmp;
clear ar_tmp;


function arSaveParOnly(ar, savepath) %#ok<INUSL>
save([savepath '/workspace_pars_only.mat'],'ar','-v7.3');
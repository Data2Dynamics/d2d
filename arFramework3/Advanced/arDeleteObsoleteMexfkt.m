% arDeleteObsoleteMexfkt
%
% This function deletes the mex-binaries in the basic folder which do not
% occur in any workspace in the results folder as ar.fkt
%
% This function is helpful, if a lot of such binaries were produced during
% model-selection by multiple changes of the def files.

function arDeleteObsoleteMexfkt

global ar

d = dir('Results');
d = d([d.isdir]==1);
folders = setdiff({d.name},{'.','..'});


fkt = cell(size(folders));
for f=1:length(folders)   
    if exist(['Results',filesep,folders{f},filesep,'workspace.mat'],'file')
        tmp = load(['Results',filesep,folders{f},filesep,'workspace.mat'],'ar');
        doload = true;
    elseif exist(['Results',filesep,folders{f},filesep,'workspace_pars_only.mat'],'file')
        doload = true;
        tmp = load(['Results',filesep,folders{f},filesep,'workspace_pars_only.mat'],'ar');
    else
        doload = false;
        fprintf('No proper workspace found in Results folder %s\n',folders{f});
    end
    if doload
        if isfield(tmp.ar,'fkt')
            fkt{f} = tmp.ar.fkt;
        else
            fkt{f} = '';
        end
    else
        fkt{f} = '';
    end        
end

fkt{end+1} = ar.fkt;

d = dir;
files = setdiff({d.name},{'.','..'});

tmp=mexext('all');
ext = {tmp.ext};
for e=1:length(ext)
    if e==1
        ind = find(~cellfun(@isempty,regexp(files,ext{e})));
    else
        ind = [ind,find(~cellfun(@isempty,regexp(files,ext{e})))];
    end
end
mexfiles = files(ind);
mexfiles = mexfiles(~cellfun(@isempty,regexp(mexfiles,'^arSimuCalcFun')));

for i=1:length(mexfiles)
    [dummy,name] = fileparts(mexfiles{i});
    if isempty(intersect(name,fkt))
        fprintf('%s is obsolete, is deleted now.\n',mexfiles{i});
        delete(mexfiles{i});
    end
end





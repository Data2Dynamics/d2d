% This function deletes the mex-binaries in the basic folder which do not
% occur in any workspace in the results folder as ar.fkt
%
% This function is helpful, if a lot of such binaries were produced during
% model-selection by multiple changes of the def files.

function arDeleteObsoleteMexfkt

d = dir('Results');
folders = setdiff({d.name},{'.','..'});


fkt = cell(size(folders));
for f=1:length(folders)
    tmp = load(['Results',filesep,folders{f},filesep,'workspace.mat'],'ar');
    fkt{f} = tmp.ar.fkt;
end

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

for i=1:length(mexfiles)
    [dummy,name] = fileparts(mexfiles{i});
    if isempty(intersect(name,fkt))
        fprintf('% is is obsolete, is deleted now.\n',mexfiles{i});
        delete(mexfiles{i});
    end
end





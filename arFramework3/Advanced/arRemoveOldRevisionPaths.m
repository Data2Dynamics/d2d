function arRemoveOldRevisionPaths

path0 = strsplit(path,pathsep);
if ispc
    old_oldrevisionPath = find(~cellfun(@isempty,regexp(path0,'Advanced\\OldRevisons')));
elseif isunix
    old_oldrevisionPath = find(~cellfun(@isempty,regexp(path0,'Advanced/OldRevisons')));
elseif ismac
    error('Please implement the corresponding command for MAC now.');
%     old_oldrevisionPath = find(~cellfun(@isempty,regexp(path0,'Advanced/OldRevisons')));
else
    error('Computer system not yet implemented');
end
    
for i=1:length(old_oldrevisionPath)
    rmpath(path0{old_oldrevisionPath(i)})
end

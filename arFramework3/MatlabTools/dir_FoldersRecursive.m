% folders = dir_FoldersRecursive(pfad)
% 
%   This functions finds all folders below the path specified by 'pfad'
% 
%   pfad      Folder name where the search should start
%             Default: pfad = pwd
% 
% Example
% folders = dir_FoldersRecursive
% folders = dir_FoldersRecursive('e:\clemens\texte')
function folders = dir_FoldersRecursive(pfad)
if(~exist('pfad','var') || isempty(pfad))
    pfad = pwd;
end

d = dir(pfad);
d = d([d.isdir]==1);
d = d(3:end); % remove '.' and '..'

fhere = strcat(pfad,filesep,{d.name});

fbelow = cell(size(fhere));
for i=1:length(fhere)    
    [fbelow{i}] = dir_FoldersRecursive(fhere{i});
end
folders = [fhere,fbelow{:}];

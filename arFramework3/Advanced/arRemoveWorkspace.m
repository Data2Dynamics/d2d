% arRemoveWorkspace
%
% Removes a workspace folder of your choice. Let's you choose from the list
% of saved workspace folders
%
% See also: arRenameWorkspace, arSave, arLoad


function arRemoveWorkspace
arCheck;

[~, filename] = fileChooser('./Results', 1, true);

remove = input(sprintf('do you really want to remove the workspace \"%s\"?\nY/[N] ', filename),'s');
if strcmpi(remove,'y')
    rmdir(['Results/' filename],'s');
    fprintf('workspace \"%s\" successfully removed!\n', filename);
end

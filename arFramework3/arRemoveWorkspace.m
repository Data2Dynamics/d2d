% Remove workspace folder of your choice
%
% Data-2-Dynamics Software
% Website: https://bitbucket.org/d2d-development/d2d-software/wiki/Home
% Contact: Andreas Raue - andreas.raue@fdm.uni-freiburg.de
% Copyright 2013 D2D Development Team. All rights reserved.

function arRemoveWorkspace
arCheck;

[~, filename] = fileChooser('./Results', 1, true);

remove = input(sprintf('do you really want to remove the workspace \"%s\"?\nY/[N] ', filename),'s');
if strcmpi(remove,'y')
    rmdir(['Results/' filename],'s');
    fprintf('workspace \"%s\" successfully removed!\n', filename);
end

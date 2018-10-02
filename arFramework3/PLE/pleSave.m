%   pleSave(ar)
% 
%   Saves the variable ar.ple to ar.config.savepath/PLE
%   in workspace result.mat
% 
%   This funciton doesn't use the info in global ar to
%   prevent that manipulations have effect on the global variable 

function pleSave(ar)
if(~exist([cd '/' ar.config.savepath '/PLE'], 'dir'))
    mkdir([cd '/' ar.config.savepath '/PLE'])
end

%warning('pleSave is deprecated. pleGlobals is now stored as ar.ple by using arSave')

ple = ar.ple;
ple.fighandel_multi = [];    % remove handle to
save([ar.config.savepath '/PLE/results.mat'], 'ple');


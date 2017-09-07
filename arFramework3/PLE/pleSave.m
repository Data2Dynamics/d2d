%  pleSave(ple)
% 
%   ple is ar.ple
% 
%   Saves the input variable ar.ple to ar.ple.savePath
%   in workspace result.mat
% 
%   This funciton doesn't use the info in global ar to
%   prevent that manipulations have effect on the global variable 
% 
% s
% 

function pleSave(ple)
if(~exist([cd '/' ple.savePath], 'dir'))
    mkdir([cd '/' ple.savePath])
end

%warning('pleSave is deprecated. pleGlobals is now stored as ar.ple by using arSave')

ple.fighandel_multi = [];    % remove handle to
save([ple.savePath '/results.mat'], 'ple');


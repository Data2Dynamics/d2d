%  pleSave(pleGlobals)
% 
%   Saves the input variable pleGlobals to pleGlobals.savePath
%   in workspace result.mat
% 
%   This funciton doesn't use the global pleGlobals to
%   prevent that manipulations have effect on the global variable 
% 
% 
% 

function    pleSave(pleGlobals)
if(~exist([cd '/' pleGlobals.savePath], 'dir'))
    mkdir([cd '/' pleGlobals.savePath])
end

pleGlobals.fighandel_multi = [];    % remove handle to
save([pleGlobals.savePath '/results.mat'], 'pleGlobals');


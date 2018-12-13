% arClearCFiles(doRemoveMex)
% 
% clear compiled c-files
%
% doRemoveMex - boolean if also clearing arSimuCalcFun* files [true]
% 
% clears complete 'Compiled' folder
%
% Example:
%   arClearCFiles

function arClearCFiles(doRemoveMex)

if(nargin==0)
    doRemoveMex = true;
end

if(doRemoveMex)
    filelist = dir;
    for j=1:length(filelist)
        indexes = strfind(filelist(j).name, 'arSimuCalcFun_');
        if(~isempty(indexes))
            delete(filelist(j).name);
        end
    end
end

% delete /Compiled folder
try %#ok<TRYNC>
    rmdir('./Compiled/', 's');
end

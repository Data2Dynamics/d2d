% clear compiled c-files

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

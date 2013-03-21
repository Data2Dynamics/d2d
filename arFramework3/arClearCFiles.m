% clear compiled c-files

function arClearCFiles

filelist = dir;

for j=1:length(filelist)
    indexes = strfind(filelist(j).name, 'arSimuCalcFun_');
    if(~isempty(indexes))
        delete(filelist(j).name);
    end
end

filelist = dir('./Compiled');

for j=1:length(filelist)
    indexes = strfind(filelist(j).name, 'mex');
    if(isempty(indexes) && ~strcmp(filelist(j).name, '.') && ~strcmp(filelist(j).name, '..'))
        delete(['./Compiled/' filelist(j).name]);
    end
end

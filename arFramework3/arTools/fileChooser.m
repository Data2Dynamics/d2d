function [out, filename] = fileChooser(filepath, default, zeigen, searchpattern)

if(~exist('default', 'var'))
    default = 1;
end
if(~exist('zeigen', 'var'))
    zeigen = true;
end
if(~exist('searchpattern', 'var'))
    searchpattern = {};
end

filelist = fileList(filepath, searchpattern);

if(length(filelist)==0)  %#ok<ISMT>
    error('No matching files found!'); 
end

out = stringListChooser(filelist, default, zeigen);
if(out ~= 0)
    filename = filelist{out};
else
    filename = '';
end


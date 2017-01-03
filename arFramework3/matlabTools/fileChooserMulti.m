function filenames = fileChooserMulti(filepath, zeigen, searchpattern)

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

out = 1;
filenames = {};
while(out~=0)
    [out, filename] = fileChooser('./Results', 0, zeigen, searchpattern);
    if(out~=0)
        filenames{end+1} = filename; %#ok<AGROW>
    end
end

function [out, filename,file_list] = fileChooser(filepath, default, zeigen, searchpattern)

if(~exist('default', 'var'))
    default = 1;
end
if(~exist('zeigen', 'var') || isempty(zeigen))
    zeigen = true;
end
if(~exist('searchpattern', 'var'))
    searchpattern = {};
end

if(zeigen==-1) % only return the file_list
    file_list = fileList(filepath, searchpattern,[],1);
    filename = [];
    out = [];
else
    
    file_list = fileList(filepath, searchpattern,[],1);
    if(length(file_list)==0)  %#ok<ISMT>
        error('No matching files found!');
    end
    
    out = stringListChooser(file_list, default, zeigen);
    if(out ~= 0)
        filename = file_list{out};
    else
        filename = '';
    end
    
end
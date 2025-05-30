% sortModus      sorting of the workspaces (passed to fileList)
%               'none' [default]
%               'chi2'
%               'checkstr'

function [out, filename,file_list] = fileChooser(filepath, default, zeigen, searchpattern,sortModus)
if(~exist('sortModus', 'var'))
    sortModus = 'none';
end

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
        error('No ''Results'' folder found! Switch your path to D2D working directory.');
    end
    
    out = stringListChooser(file_list, default, zeigen, sortModus);
    if(out ~= 0)
        filename = file_list{out};
    else
        filename = '';
    end
    
end
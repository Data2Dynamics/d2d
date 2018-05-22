%   This function can be used to list all filenames recursively from the
%   working directory (incl. the path)

function all_files = list_files_recursive

check_folders = {strcat(pwd)};
all_files = cell(0);

while ~isempty(check_folders)
    tmpfolder = check_folders{1};
%     disp(tmpfolder)
    [files,folders] = list_core(tmpfolder);
    all_files = union(all_files,files);
    check_folders = setdiff(union(check_folders,folders),tmpfolder);    
end





function [files,folders] = list_core(ordner)

d = dir(ordner);
folders = {d.name};

folders = folders(find([d.isdir]));
folders = setdiff(folders,{'.','..'});
folders = strcat(strcat(ordner,filesep),folders);


files = {d.name};
files = files(find(~[ d.isdir]));
files = setdiff(files,{'.','..'});
files = strcat(strcat(ordner,filesep),files);

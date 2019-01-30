%   [all_files,all_folders] = list_files_recursive(max_depth,mode)
% 
%   This function can be used to list all filenames recursively from the
%   working directory (incl. the path)
% 
% 
%   mode            can be used to speed up:
%                   if 1 or 'files', then only files will be searched
%                   if 2 or 'folders' then only folders will be searched
%                   default: 0 [this means both]
% 
%   all_files       all filenames including path
%   all_folders     all foldernames including path
% 
% 
%   max_depth       the maximal depth, the current directory has depth 0
%                   default: Inf
% 
% 
% Examples:
% 
% list_files_recursive 
%       lists all files recursively
% 
% all_files = list_files_recursive(0)
%       lists all files in the current folder
%  
% all_files = list_files_recursive(1)
%       lists all files in the current folder and in all subfolders
% 
% [~,f]=list_files_recursive(1) % all folders in pwd and one subfolder deeper


function [all_files,all_folders] = list_files_recursive(max_depth,mode)
if ~exist('max_depth','var') || isempty(max_depth)
    max_depth = Inf;
end
if ~exist('mode','var') || isempty(mode)
    mode = 0;
elseif strcmpi(mode,'files') % has no effect in the current implementation
    mode = 1;
elseif strcmpi(mode,'folders')
    mode = 2;
else
    mode
    error('mode invalid')
end


check_folders = {strcat(pwd)};
all_files = cell(0);
all_folders = cell(0);

depth_pwd = length(strfind(pwd,filesep)); % number of filesep (\ or /)


while ~isempty(check_folders)
    tmpfolder = check_folders{1};
    %     disp(tmpfolder)
    [files,folders] = list_core(tmpfolder);
    depth = length(strfind(tmpfolder,filesep)) - depth_pwd;

    if depth <= max_depth
        if mode~=2
            all_files = union(all_files,files);
        end
        all_folders = union(all_folders,folders);
            
        check_folders = setdiff(union(check_folders,folders),tmpfolder);
    else
        check_folders = setdiff(check_folders,tmpfolder);        
    end
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

% [has_git, is_repo] = arCheckGit(repo_path)
% 
% Check if git is available on system and if d2d is installed as a valid git repository
%
%   repo_path   the path where D2D is located as e.g. stored in ar.info.ar_path
%   has_git     logical indicating whether git is available
%   is_repo     logical indicating whether D2D is downloaded via git
% 
% The shell command "git status" is called.
% 
% See also arCheckVersion, arUpdateD2D

function [has_git, is_repo] = arCheckGit(repo_path)

if(~exist('repo_path','var'))
	repo_path = pwd;
end

% check if git exists on system and suppress output
if(isunix)
    has_git = system('which git >/dev/null 2>&1')==0;
else
    has_git = system('where git >nul 2>&1')==0;
end

% check if code is organized as a git reopository (false for zip/tar.gz downloads)
if(has_git)
	old_dir = pwd;
	cd(repo_path)
    try
        if(isunix)
            is_repo = system('git status >/dev/null 2>&1')==0;
        else
            is_repo = system('git status >nul 2>&1')==0;
        end
        cd(old_dir)
    catch err
        cd(old_dir)
        rethrow(err)
    end
else
	is_repo = false;
end



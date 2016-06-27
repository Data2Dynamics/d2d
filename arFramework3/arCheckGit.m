% Check if git is available on system and if d2d is installed as a valid git repository
%
% [has_git, is_repo] = arCheckGit(repo_path)

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
    if(isunix)
        is_repo = system('git status >/dev/null 2>&1')==0;
    else
        is_repo = system('git status >nul 2>&1')==0;
    end
    cd(old_dir)
else
	is_repo = false;
end



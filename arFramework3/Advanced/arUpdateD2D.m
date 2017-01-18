% update D2D to current revision on github

function arUpdateD2D

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

ar_path = ar.info.ar_path;
[has_git, is_repo] = arCheckGit(ar_path);

flag = -1;
if arCheckVersion(true)~=1
    if(isunix)
        library_path = getenv('LD_LIBRARY_PATH');
        setenv('LD_LIBRARY_PATH', '');
    end
    if (has_git && is_repo)
        old_path = pwd;
        
        cd(ar_path)
        flag = system('git pull origin master');
        cd(old_path)
    end
    if(isunix)
        setenv('LD_LIBRARY_PATH', library_path);
    end
else
    flag = 2;
end

switch flag
    case 0
        arFprintf(2,'D2D was updated successfully\n')
    case 2
        arFprintf(2,'D2D is already up-to-date\n')
    otherwise
        arFprintf(2,'There was a problem while updating D2D\n')
end

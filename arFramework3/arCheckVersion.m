% check D2D version and compare with current revision on github

function arCheckVersion

global ar  

[has_git, is_repo] = arCheckGit(ar.info.ar_path);

if(has_git && is_repo)
    % find current revision of d2d
    old_path = pwd;
    cd(ar.info.ar_path)
    [~, cmdout] = system('git rev-parse HEAD');
    cd(old_path)
    clear old_path;

    ar.info.revision = deblank(cmdout);
    
    % get current SHA from github and compare with installed revision
    try
        if ( exist('webread', 'file')==2 )
            gh_data = webread('https://api.github.com/repos/Data2Dynamics/d2d/git/refs/heads/master');
            if(~isempty(gh_data) && ~strcmp(deblank(gh_data.object.sha),ar.info.revision))
                warning( 'There is a newer version available on github! Please check http://www.data2dynamics.org for updates.' );
            end
        else
            if ( exist('urlread', 'file')==2 )
                % Code path for older versions of MATLAB
                try
                    gh_data = urlread('https://api.github.com/repos/Data2Dynamics/d2d/git/refs/heads/master');
                    gh_data = gh_data(strfind(gh_data, '"sha":')+7:end);
                    br      = strfind(gh_data, '"');
                    gh_data = gh_data(1:br(1)-1);
                    clear br;
                    if(~isempty(gh_data) && ~strcmp(deblank(gh_data),ar.info.revision))
                        warning( 'There is a newer version available on github! Please check http://www.data2dynamics.org for updates.' );
                    end
                catch
                end 
            end
        end
    catch err
        warning('Could not connect to github. Version checking cancelled.')
    end
else
    ar.info.revision = [];
end

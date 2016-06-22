% Initialize and clear workspace of framework
%
% Data-2-Dynamics Software
% Website: https://bitbucket.org/d2d-development/d2d-software/wiki/Home
% Contact: Andreas Raue - andreas.raue@fdm.uni-freiburg.de
% Copyright 2013 D2D Development Team. All rights reserved.

if(~arCheck)
    return;
end

global ar

warning('off', 'symbolic:sym:sym:DeprecateExpressions')

ar = struct([]);
ar(1).stop = 0;
ar.fevals = 0; 

arInitUser;

ar.info.initTime = now;
[ar.info.def_version_code, ar.info.c_version_code] = arGetVersion;
arFprintf(1, 'Data 2 Dynamics Software\n');
arFprintf(1, '(arFramework3, def-version %i, c-version %s)\n', ...
    ar.info.def_version_code, ar.info.c_version_code);
arFprintf(1, 'Website: http://www.data2dynamics.org\n');
arFprintf(1, 'Contact: Andreas Raue - andreas.raue@fdm.uni-freiburg.de\n');
arFprintf(1, 'Copyright 2015 D2D Development Team. All rights reserved.\n\n');

ar.checksum = [];

if(isunix)
    % check if git exists on system
    has_git = system('which git > /dev/null')==0;
else
    has_git = false;
end

if(has_git)
    % find current revision of d2d
    ar_path = fileparts(which('arInit.m'));
    ar.info.ar_path = ar_path;
    old_path = pwd;
    cd(ar_path)
    [~,cmdout] = system('git rev-parse HEAD');
    cd(old_path)

    ar.info.revision = deblank(cmdout);
    cmdout = [];
    
    % get current SHA from github and compare with installed revision
    try
        if ( exist('webread') ) %#ok
            gh_data = webread('https://api.github.com/repos/Data2Dynamics/d2d/git/refs/heads/master');
            if(~isempty(gh_data) && ~strcmp(deblank(gh_data.object.sha),ar.info.revision))
                warning( 'There is a newer version available on github! Please check http://www.data2dynamics.org for updates.' );
            end
        else
            if ( exist('urlread') ) %#ok
                % Code path for older versions of MATLAB
                try
                    gh_data = urlread('https://api.github.com/repos/Data2Dynamics/d2d/git/refs/heads/master');
                    gh_data = gh_data(strfind(gh_data, '"sha":')+7:end);
                    br      = strfind(gh_data, '"');
                    gh_data = gh_data(1:br(1)-1);
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

ar = arInitFields(ar);

% check licenses
if ( ~license('test', 'Symbolic_Toolbox') )
    warning( 'D2D requires a license for the MathWorks symbolic math toolbox. It is unlikely that D2D will work.' );
end
if ( ~license('test', 'Optimization_Toolbox') )
    warning( 'No license found for optimization toolbox. If fitting is required, obtain a license or switch optimization method (e.g. ar.config.optimizer=3).' );
end

ar = orderfields(ar);
ar.info = orderfields(ar.info);
ar.config = orderfields(ar.config);
ar.ppl = orderfields(ar.ppl);

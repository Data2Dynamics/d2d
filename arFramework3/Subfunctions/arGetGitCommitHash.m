% Get the hash of the git commit for saving in ar.info.gitCommitHash
%
% The D2d version used to build the ar struct can be obtained by
% https://github.com/Data2Dynamics/d2d/tree/hash
% where 'hash' has to be replaced by the string saved in ar.info.gitCommitHash
%
% commit_hash = arGetGitVersion(repo_path)

function commit_hash = arGetGitCommitHash(repo_path)

if(~exist('repo_path','var'))
    repo_path = fileparts(which('arInit.m'));
end

old_dir = pwd;
cd(repo_path)
try
    if(isunix)
        [~, commit_hash] = system('git rev-parse HEAD 2>/dev/null');
    else
        [~, commit_hash] = system('git rev-parse HEAD 2>nul');
    end
    commit_hash(commit_hash<30) = [];  % removes line breaks 
catch ME
    commit_hash = '';
    warning('Error in fetching git commit hash. No git version set. Error meassage:')
    fprintf(getReport(ME));
end

cd(old_dir)
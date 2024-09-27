% arSetParallelThreads([n], [silent])
%
% Set how many parallel threads of the machine should be used
%
% n         number of threads (default = 2*number of cores of machine)
% silent    if 'true' supresses output (default = 'false', output is printed)
%
% Example:
%   arSetParallelThreads
%   >> requesting 8 thread(s) for 23 task(s) on 4 core(s).
%
%   arSetParallelThreads(1)
%   >> requesting 1 thread(s) for 23 task(s) on 4 core(s).


function arSetParallelThreads(n, silent)

global ar

% get correct number of cores and max threads
ar.config.nCore = feature('numCores');
ar.config.nMaxThreads = 2*ar.config.nCore;

% set default values
if ~exist('silent','var') || isempty(silent)
    silent = false;
end
if ~exist('n','var') || isempty(n)
    n = ar.config.nMaxThreads;
end

if ar.config.useParallel == 0
    if ~silent
        fprintf('parallel computing is disabled, reset requested threads to 1.\n');
        fprintf('enable parallel computing with "ar.config.useParallel=1".\n');
    end
    n = 1;

else
    if ~silent
        fprintf('requesting %i thread(s) for %i task(s) on %i core(s).\n', ...
            n, ar.config.nTasks, ar.config.nCore);
    end
    if n > ar.config.nMaxThreads
        if ~silent
            fprintf('request exceeds thread limit, reset requested threads to %i.\n', ...
                ar.config.nMaxThreads);
        end
        n = ar.config.nMaxThreads;
    end
    if n > ar.config.nTasks
        if~silent
            fprintf('request exceeds number of tasks, reset requested threads to %i.\n', ...
                ar.config.nTasks);
        end
        n = ar.config.nTasks;
    end
end

ar.config.nParallel = n;
arLink(true);

end
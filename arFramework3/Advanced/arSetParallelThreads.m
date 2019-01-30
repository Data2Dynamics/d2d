% arSetParallelThreads([n],[silent])
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

ar.config.nCore = feature('numCores');

if(~exist('n','var'))
    n = 2*ar.config.nCore;
end
if(~exist('silent','var'))
    silent = false;
end

if(~silent)
    fprintf(1,'requesting %i thread(s) for %i task(s) on %i core(s).\n', ...
        n, ar.config.nTasks, ar.config.nCore);
end

if(n>ar.config.nTasks)
    if(~silent)
        fprintf(1,'less tasks than %i cores, reset requested threads to %i.\n', ...
            ar.config.nCore, ar.config.nTasks);
    end
    n = ar.config.nTasks;
end

ar.config.nParallel = n;
arLink(true);

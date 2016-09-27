% set how many parallel thread n should be executed
%
% arSetParallelThreads(n)
%
% n default = number of cores of machine

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

if(ar.config.nParallel>ar.config.nTasks)
    if(~silent)
        fprintf(1,'less tasks than %i cores, reset requested threads to %i.\n', ...
            ar.config.nCore, ar.config.nTasks);
    end
    n = ar.config.nTasks;
end

ar.config.nParallel = n;
arLink(true);

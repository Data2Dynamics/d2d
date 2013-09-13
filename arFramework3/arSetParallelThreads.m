% set how many parallel thread n should be executed
%
% arSetParallelThreads(n)
%
% n default = number of cores of machine

function arSetParallelThreads(n)

global ar

ar.config.nCore = feature('numCores');

if(~exist('n','var'))
    n = ar.config.nCore;
end

fprintf('requesting %i threads for %i tasks on %i core machine.\n', ...
    n, ar.config.nTasks, ar.config.nCore);

ar.config.nParallel = n;
arLink(true);

if(ar.config.nParallel>ar.config.nTasks)
    fprintf('less tasks than %i cores, reset requested threads to %i.\n', ...
        ar.config.nTasks, ar.config.nCore, ar.config.nTasks);
end

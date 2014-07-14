% arCheckParallelSpeedUp(nmax, nrep, sensi, randomseed)

function arParallelSpeedUpTest(nmax, nrep, sensi, randomseed)

global ar

if(~exist('nmax','var') || isempty(nmax))
    nmax = ar.config.nCore + ceil(ar.config.nCore/2);
    if(nmax > ar.config.nTasks)
        nmax = ar.config.nTasks;
    end
end
if(~exist('nrep','var') || isempty(nrep))
    nrep = 100;
end
if(~exist('sensi','var') || isempty(sensi))
    sensi = false;
end
if(~exist('randomseed','var'))
    randomseed = 1;
end

nreset = ar.config.nParallel;

% make vector
if(length(nmax)==1)
    nmax = 1:nmax;
end

% make sure 1 is included, an nmax is sorted
nmax = unique([1 nmax]);
T = nan(nrep,length(nmax));
neffective = nan(1,length(nmax));
for jn = 1:length(nmax)
    arSetParallelThreads(nmax(jn));
    
    neffective(jn) = 0;
    for jt = 1:length(ar.config.threads);
        neffective(jn) = max([neffective(jn) ar.config.threads(jt).n]);
    end
    
    if(ar.config.nParallel>ar.config.nTasks)
        break;
    end
    
    arChi2LHS(nrep, sensi, randomseed, true)
    T(:,jn) = ar.timing;
end
ar.config.nParallel = nreset;

nTasks = ar.config.nTasks;
nCore = ar.config.nCore;

savestr = ['arCheckParallelSpeedUp_' datestr(now, 30) '.mat'];
save(savestr, 'T', 'nTasks', 'nCore', 'nmax', 'neffective');

%% plot

Ttmp = bsxfun(@rdivide, T(:,1), T(:,2:end));

nthreads = nmax(2:end);
tspeedup = nTasks./neffective(2:end);

figure(1); clf;
boxplot(Ttmp, 'orientation', 'horizontal', ...
    'plotstyle', 'compact', 'colors', 'k');
hold on
h1 = plot(tspeedup, 1:length(nthreads), 'r*-');
xlim([1 max([Ttmp(:); tspeedup(:)])*1.1])
h2 = plot(xlim, [1 1]*((nCore-1)+0.5), 'k--');
hold off
ylabel('number of threads');
xlabel('fold speed up');
title('speed up of parallel computation compared to sequential computation');
set(gca, 'yTick', 1:length(nthreads));
set(gca, 'yTickLabel', nthreads);
set(gca, 'XScale', 'log');
legend([h1 h2],'expected speed up','machine core limit','Location','SouthEast')


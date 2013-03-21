function arTestParallelSpeedUp(N, sensis)

global ar

ar.config.useParallel = false;
arChi2LHS(N, sensis, true);
ar.T1 = ar.timing;

ar.config.useParallel = true;
arChi2s(ar.ps, sensis, true);
ar.T2 = ar.timing;

%% plot

figure(2)
hist(log2(ar.T1./ar.T2), 20);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.7 .7 .7]);
set(h,'EdgeColor',[0 0 0]);


% title(sprintf('parallelization acceleration 1-core vs. %i-cores',maxNumCompThreads))
title({'parallelization acceleration',sprintf('1-core vs. %i-cores',feature('numCores'))})
xlabel('log-2 fold reduction in computation time')
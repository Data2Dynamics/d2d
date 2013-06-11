function arPlotChi2s

global ar

figure(1)

[chi2s, isorted] = sort(ar.chi2s);
exitflag = ar.exitflag(isorted);

chi2min = min([chi2s ar.chi2fit]);
ar.chi2s_sorted = chi2s;
ar.ps_sorted = ar.ps(isorted,:);

semilogy(chi2s - chi2min + 1, '--');
hold on
h = semilogy(find(exitflag>0), chi2s(exitflag>0) - chi2min + 1, 'o');
plot(xlim, [1 1], 'k--');
hold off
if(ar.config.fiterrors == 1)
    ylabel('-2*log(L) + const');
else
    ylabel('\chi^2');
end
xlabel('evaluations sorted');
title(sprintf('%i evaluations in total, %i without errors, %i converged', ...
    length(exitflag), sum(~isnan(chi2s)) ,sum(exitflag>0)));
legend(h, 'converged evaluations');

figure(2)
hist(ar.timing, 50);
xlabel('time / sec.');
title(sprintf('total time for %i evaluations %s', ...
    length(ar.chi2s), secToHMS(sum(ar.timing(~isnan(ar.timing))))));

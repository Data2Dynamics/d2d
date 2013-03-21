function arPlotChi2s

global ar

figure(1)

[chi2s, isorted] = sort(ar.chi2s);
exitflag = ar.exitflag(isorted);

ar.chi2s_sorted = chi2s;
ar.ps_sorted = ar.ps(isorted,:);

% bounds = sum(bsxfun(@gt, ar.ps_sorted, ar.ub-0.1) | ...
%     bsxfun(@lt, ar.ps_sorted, ar.lb+0.1),2)';

if(ar.config.fiterrors == 1)    
%     plot(sort(2*ar.ndata*log(sqrt(2*pi)) + ar.chi2s), 'o--');
    semilogy(chi2s - ar.chi2fit + 1, '--');
    hold on
    h = semilogy(find(exitflag>0), chi2s(exitflag>0) - ar.chi2fit + 1, 'o');
%     h2 = semilogy(find(bounds>0), chi2s(bounds>0) - ar.chi2fit + 1, 'r*');
    plot(xlim, [1 1], 'k--');
    hold off
    ylabel('-2*log(L) + const');
else
    semilogy(chi2s, '--');
    hold on
    h = semilogy(find(exitflag>0), chi2s(exitflag>0), 'o');
    plot(xlim, [ar.chi2fit ar.chi2fit], 'k--');
    hold off
    ylabel('\chi^2');
end
xlabel('fits sorted');
title(sprintf('%i operations in total, %i without errors, %i converged', ...
    length(exitflag), sum(~isnan(chi2s)) ,sum(exitflag>0)));
legend(h, 'converged operations');
% legend([h h2], 'converged fits', 'at bounds');

% T = clusterdata(log10(chi2s(~isnan(chi2s))'), 1)
% Y = pdist(log10(chi2s(~isnan(chi2s))'));
% Z = linkage(Y);
% figure(3)
% [H,T] = dendrogram(Z,'colorthreshold','default');

figure(2)
hist(ar.timing, 50);
xlabel('fit time / sec.');
title(sprintf('total time for %i operations %s', ...
    length(ar.chi2s), secToHMS(sum(ar.timing(~isnan(ar.timing))))));

if(~isempty(ar.fun_evals))
    figure(3)
    hist(ar.fun_evals, 50);
end
xlabel('number of function evaluations');
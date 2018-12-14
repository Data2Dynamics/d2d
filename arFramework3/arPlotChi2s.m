% arPlotChi2s
%
% This function creates one figure and plots information about chi2s and 
% runtime of the multistart. It is intended to be used after arChi2LHS.
% 
% see also arPlotFits arChi2LHS arFitLHS

function arPlotChi2s

global ar

if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

figure(1)

subplot(3,1,[1 2]);

[chi2s, isorted] = sort(ar.chi2s);
ar.chi2s_sorted = chi2s;
ar.chi2sconstr_sorted = ar.chi2sconstr(isorted);
ar.ps_sorted = ar.ps(isorted,:);
exitflag = ar.exitflag(isorted);

chi2min = min([chi2s arGetMerit('chi2fit')]);
chi2minconstr = min([ar.chi2sconstr_sorted arGetMerit('chi2constr')]);

h = [];
semilogy(chi2s - chi2min + 1, '--');
hold on
h(2) = semilogy(ar.chi2sconstr_sorted - chi2minconstr + 1, 'r--');
h(1) = semilogy(find(exitflag>0), chi2s(exitflag>0) - chi2min + 1, 'o');
plot(xlim, [1 1], 'k--');
hold off
if ar.config.fiterrors==1 || (ar.config.fiterrors==0 && sum(ar.qFit==1 & ar.qError==1)>0)
    ylabel('-2 log(L) + const');
else
    ylabel('\chi^2');
end
xlabel('run index (sorted by likelihood)');
title(sprintf('%i evaluations in total, %i without errors, %i converged', ...
    length(exitflag), sum(~isnan(chi2s)) ,sum(exitflag>0)));
legend(h, {'converged evaluations', 'constraint violation'});

subplot(3,1,3);
hist(ar.timing, 50);
xlabel('time / sec.');
title(sprintf('total time for %i evaluations %s', ...
    length(ar.chi2s), secToHMS(sum(ar.timing(~isnan(ar.timing))))));

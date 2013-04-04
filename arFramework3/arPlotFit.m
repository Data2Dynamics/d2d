% Plot fit of model parameters to data using Levenberg-Marquardt
%
% arPlotFitLM

function arPlotFit

global ar

figure(2)

chi2s = ar.fit.chi2_hist;
xs = 1:sum(~isnan(chi2s));

subplot(3,2,1)
plot(xs, chi2s(xs), '-')
title('likelihood improvement')

subplot(3,2,3)
plot(xs(1:end-1), log10(abs(-diff(chi2s(xs)))), '-')
hold on
plot(xlim, log10([ar.config.optim.TolFun ar.config.optim.TolFun]), 'r--');
hold off
title('log10 relativ likelihood improvement')

subplot(3,2,5)
plot(xs, log10(ar.fit.maxstepsize_hist(xs)), 'r-')
hold on
plot(xs, log10(ar.fit.stepsize_hist(xs)), '-')
plot(xlim, log10([ar.config.optim.TolX ar.config.optim.TolX]), 'r--');
hold off
title('log10 stepsize')
xlabel('iterations')

subplot(3,2,2)
plot(xs, ar.fit.p_hist(xs,:))
title('parameters')

subplot(3,2,4)
plot(xs, bsxfun(@minus,ar.fit.p_hist(1,:),ar.fit.p_hist(xs,:)))
title('parameter changes relative to start')

subplot(3,2,6)
plot(xs(1:end-1), diff(ar.fit.p_hist(xs,:))) 
title('relative parameter changes')
xlabel('iterations')
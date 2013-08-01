% Plot fit

function arPlotFit(index, qp)

global ar

if(~exist('index','var'))
    fit = ar.fit;
else
    fit = ar.fit_hist(index).hist;
end

if(~exist('qp','var'))
    qp = true(size(ar.p));
end

figure(2); clf;

chi2s = fit.chi2_hist;
opti = fit.opti_hist;
xs = 1:sum(~isnan(chi2s));

subplot(4,2,1)
plot(xs-1, log10(chi2s(xs) - min(chi2s(xs)) + 1), '-')
title('log10 likelihood improvement')
xlim([0 length(xs)-1])

subplot(4,2,3)
plot(xs-1, log10(fit.stepsize_hist(xs)), '-')
hold on
plot(xs-1, log10(fit.maxstepsize_hist(xs)), 'r--')
plot([0 length(xs)], log10([ar.config.optim.TolX ar.config.optim.TolX]), 'r--');
hold off
title('log10 stepsize')
xlabel('iterations')
xlim([0 length(xs)-1])

subplot(4,2,5)
plot(xs(2:end)-1, log10(abs(-diff(chi2s(xs)))), '-')
hold on
plot(xlim, log10([ar.config.optim.TolFun ar.config.optim.TolFun]), 'r--');
hold off
title('log10 relativ likelihood improvement')
xlim([0 length(xs)-1])

subplot(4,2,7)
plot(xs-1, log10(opti(xs)), '-')
hold on
plot([0 length(xs)], [-6 -6], 'r--');
hold off
title('log10 first order optimality criterion')
xlabel('iterations')
xlim([0 length(xs)-1])

subplot(4,2,2)
plot(xs-1, fit.p_hist(xs,qp))
title('parameters - absolute value')
xlim([0 length(xs)-1])

subplot(4,2,4)
ps = fit.p_hist(xs,qp);
lb = ar.lb(qp);
ub = ar.ub(qp);
ps = bsxfun(@rdivide, bsxfun(@minus, ps, lb), ub-lb);
plot(xs-1, ps)
title('parameters - relative to normalized bounds')
xlim([0 length(xs)-1])

subplot(4,2,6)
plot(xs-1, bsxfun(@minus,fit.p_hist(1,qp),fit.p_hist(xs,qp)))
title('parameters - changes relative to fit start')
xlim([0 length(xs)-1])

subplot(4,2,8)
plot(xs(2:end)-1, diff(fit.p_hist(xs,qp))) 
title('parameters - differential changes')
xlabel('iterations')
xlim([0 length(xs)-1])
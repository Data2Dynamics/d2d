% Plot convergence of a single fit

function arPlotFit(index, qp)

global ar

if(~exist('index','var'))
    fit = ar.fit;
    fprintf('%s\n','current');
else
    fit = ar.fit_hist(index).hist;
    fprintf('%s\n',ar.fit_hist(index).name);
end

if(~exist('qp','var'))
    qp = ar.qFit==1;
end

chi2s = fit.chi2_hist;
if(isfield(fit,'constr_hist'))
    constrs = fit.constr_hist;
    chi2sconstrs = chi2s + constrs;
end
opti = fit.opti_hist;
xs = 1:sum(~isnan(chi2s));
if(length(xs)<2)
    fprintf('no successfull iterations\n');
    return;
end

figure(2); clf;

subplot(4,2,1)
h = [];
labels = {};
if(ar.ndata>0)
    h(end+1) = plot(xs-1, log10(chi2s(xs) - min(chi2s(xs)) + 1), 'k-');
    labels{end+1} = 'likelihood';
    hold on
end
if(ar.nconstr>0)
    h(end+1) = plot(xs-1, log10(constrs(xs) - min(constrs(xs)) + 1), 'b-');
    labels{end+1} = 'constraints';
    hold on
end
if(ar.ndata>0 && ar.nconstr>0)
    h(end+1) = plot(xs-1, log10(chi2sconstrs(xs) - min(chi2sconstrs(xs)) + 1), 'm-');
    labels{end+1} = 'likelihood + constraints';
end
hold off
legend(h, labels);
title('log10 improvement')
xlim([0 length(xs)-1])

subplot(4,2,3)
if(ar.ndata>0)
    diffchi2s = log10(-diff(chi2s(xs)));
    qreal = imag(diffchi2s)==0;
    diffchi2s(~qreal) = nan;
    plot(xs(2:end)-1, diffchi2s, 'k-');
    hold on
end
if(ar.nconstr>0)
    diffconstr = log10(-diff(constrs(xs)));
    qreal = imag(diffconstr)==0;
    diffconstr(~qreal) = nan;
    plot(xs(2:end)-1, diffconstr, 'b-');
    hold on
end
if(ar.ndata>0 && ar.nconstr>0)
    diffchi2sconstrs = log10(-diff(chi2sconstrs(xs)));
    qreal = imag(diffchi2sconstrs)==0;
    diffchi2sconstrs(~qreal) = nan;
    plot(xs(2:end)-1, diffchi2sconstrs, 'm-');
end
plot(xlim, log10([ar.config.optim.TolFun ar.config.optim.TolFun]), 'r--');
hold off

title('log10 relativ improvement')
xlim([0 length(xs)-1])

subplot(4,2,5)
plot(xs-1, log10(opti(xs)), '-')
hold on
plot([0 length(xs)], [-6 -6], 'r--');
hold off
title('log10 first order optimality criterion')
xlabel('iterations')
xlim([0 length(xs)-1])

subplot(4,2,7)
plot(xs-1, log10(fit.stepsize_hist(xs)), '-')
hold on
plot(xs-1, log10(fit.maxstepsize_hist(xs)), 'r--')
plot([0 length(xs)], log10([ar.config.optim.TolX ar.config.optim.TolX]), 'r--');
hold off
title('log10 stepsize')
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

return

figure(3); clf;
% xs = 100:500;
% xs = 100:200;
% xs = 200:400;
% xs = 1:100;
% xs = 100:300;
% xs = 250:450;
ps = fit.p_hist(xs,:);
sumonbound = sum(bsxfun(@eq, ps(2:end,:), ar.lb) | bsxfun(@eq, ps(2:end,:), ar.ub),2);
psdiff = diff(ps);
xs2 = xs(2:end)-1;
[nrows, ncols] = arNtoColsAndRows(sum(qp)+1);
subcount = 1;
for j=find(qp)
    subplot(nrows, ncols,subcount)
    plot(xs2, psdiff(:,j));
    hold on
    qonbound = ps(2:end,j)==ar.lb(j) | ps(2:end,j)==ar.ub(j);
    if(sum(qonbound)>0)
        plot(xs2(qonbound), psdiff(qonbound,j), 'rx');
    end
    hold off
    title(strrep(ar.pLabel{j}, '_', '\_'));
	arSpacedAxisLimits
    subcount = subcount + 1;
end
subplot(nrows, ncols,subcount)
plot(xs2, sumonbound);

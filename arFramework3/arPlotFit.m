% Plot fit of model parameters to data using Levenberg-Marquardt
%
% arPlotFitLM

function arPlotFit

global ar

figure(2)

qnonnan = ~isnan(ar.fit.chi2_hist);

chi2s = ar.fit.chi2_hist;

subplot(4,1,1)
plot(0:(sum(qnonnan)-1), chi2s(qnonnan), '-')
xlim([0 sum(qnonnan)-1])
title('fit improvement')

subplot(4,1,2)
% plot(0:(sum(qnonnan)-1), chi2s(qnonnan), '-')
semilogy(0:(sum(qnonnan)-2), -diff(chi2s(qnonnan)), '-')
xlim([0 sum(qnonnan)-1])
title('fit improvement')

subplot(4,1,3)
semilogy(0:(sum(qnonnan)-1), ar.fit.lambda_hist(qnonnan), '-')
xlim([0 sum(qnonnan)-1])
title('lambda')

subplot(4,1,4)
if(isfield(ar,'pTrue'))
%     plot(0:(sum(qnonnan)-1), ones(sum(qnonnan),1)*ar.pTrue-ar.fit.p_hist(1:sum(qnonnan),:))
    plot(0:(sum(qnonnan)-1), ar.fit.p_hist(1:sum(qnonnan),:))
else
    plot(0:(sum(qnonnan)-1), ones(sum(qnonnan),1)*ar.fit.p_hist(1,:)-ar.fit.p_hist(1:sum(qnonnan),:))
end
xlim([0 sum(qnonnan)-1])
title('parameter changes')
xlabel('iterations')

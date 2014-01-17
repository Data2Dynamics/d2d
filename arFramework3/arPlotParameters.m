

function arPlotParameters(nbest, qsub)

global ar

if(~exist('nbest','var') || isempty(nbest))
    nbest = size(ar.ps_sorted,1);
end
if(~exist('qsub','var'))
    qsub = true(size(ar.p));
end

[chi2s, isorted] = sort(ar.chi2s);
ar.chi2s_sorted = chi2s;
ar.ps_sorted = ar.ps(isorted,:);

ps = ar.ps_sorted(:,qsub==1);

figure(1)

colors = jet(nbest);
for jp = nbest:-1:1
    plot(ps(jp,:), 1:sum(qsub), 'Color', colors(jp,:))
    hold on
end
hold off

ylim([0 length(ar.p(qsub==1))+1])
set(gca,'YTick',1:length(ar.p(qsub==1)));
% set(gca,'YTickLabel',strrep(ar.pLabel(qsub==1),'_','\_'));
set(gca,'YTickLabel',ar.pLabel(qsub==1));
% xticklabel_rotate

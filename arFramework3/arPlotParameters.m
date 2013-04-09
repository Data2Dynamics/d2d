

function arPlotParameters(nbest)

global ar

if(~exist('nbest','var'))
    nbest = size(ar.ps_sorted,1);
end

[chi2s, isorted] = sort(ar.chi2s);
ar.chi2s_sorted = chi2s;
ar.ps_sorted = ar.ps(isorted,:);

ps = ar.ps_sorted;

figure(1)

colors = jet(nbest);
for jp = nbest:-1:1
    plot(ps(jp,:), 'Color', colors(jp,:))
    hold on
end
hold off

xlim([0 length(ar.p)+1])
set(gca,'XTick',1:length(ar.p));
set(gca,'XTickLabel',strrep(ar.pLabel,'_','\_'));
xticklabel_rotate

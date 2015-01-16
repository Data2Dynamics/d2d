% arPlotParameters(nbest, jks)

function arPlotParameters(nbest, jks)

global ar

if(~exist('nbest','var') || isempty(nbest))
    nbest = size(ar.ps_sorted,1);
end
if(~exist('jks','var'))
    jks = 1:length(ar.p);
end

[chi2s, isorted] = sort(ar.chi2s);
ar.chi2s_sorted = chi2s;
ar.ps_sorted = ar.ps(isorted,:);

ps = ar.ps_sorted(:,jks);

figure(1)

colors = jet(nbest);
for jp = nbest:-1:1
    plot(ps(jp,:), 1:length(jks), 'Color', colors(jp,:))
    hold on
end
hold off

ylim([0 length(ar.p(jks))+1])
set(gca,'YTick',1:length(ar.p(jks)));
set(gca,'YTickLabel', arNameTrafo(ar.pLabel(jks)));


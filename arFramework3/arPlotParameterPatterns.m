function arPlotParameterPatterns(ps, jks)

global ar

if(~exist('nbest','var'))
    nbest = size(ps,1);
end
if(~exist('jks','var'))
    jks = 1:length(ar.p);
end

figure(1)
C = jet(nbest);
for j=1:nbest
    plot(ps(j,jks), 1:length(jks), 'Color', C(j,:));
    hold on
end
hold off
title(sprintf('parameter changes for %i best fits', nbest));
xlabel('parameter values')
set(gca, 'YTick', 1:length(jks));
set(gca, 'YTickLabel', arNameTrafo(ar.pLabel(jks)));
ylim([0 length(jks)+1])
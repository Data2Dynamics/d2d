function arCompareFit(indexes)

global ar

if(~exist('indexes','var'))
    indexes = 1:length(ar.fit_hist);
end

minchi2 = Inf;
for j=1:length(indexes)
    minchi2 = min([minchi2 ar.fit_hist(j).hist.chi2_hist]);
end

figure(1)
h = nan(1,length(indexes));
labels = {};
for j=1:length(indexes)
    C = arLineMarkersAndColors(j,[],[],'-');
    h(j) = semilogy(ar.fit_hist(j).hist.chi2_hist + 1 - minchi2, C{:});
    labels{j} = ar.fit_hist(j).name; %#ok<AGROW>
    hold on
end
hold off

legend(h, strrep(labels, '_', '\_'));
xlabel('fit iteration');
ylabel('likelihood');
grid on
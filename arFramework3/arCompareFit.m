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
colors = jet(length(indexes));
colors = bsxfun(@rdivide, colors, sqrt(sum(colors.^2,2)));
labels = {};
for j=1:length(indexes)
    h(j) = semilogy(ar.fit_hist(j).hist.chi2_hist + 1 - minchi2, 'o-', 'Color', colors(j,:), ...
        'MarkerFaceColor','w', ...
        'LineWidth',1, 'MarkerSize',4);
    labels{j} = ar.fit_hist(j).name; %#ok<AGROW>
    hold on
end
hold off

legend(h, strrep(labels, '_', '\_'));
xlabel('fit iteration');
ylabel('likelihood');

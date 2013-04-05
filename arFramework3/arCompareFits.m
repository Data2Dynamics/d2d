function arCompareFits(filenames)

if(nargin==0)
    filenames = fileChooserMulti('./Results', true);
end

minchi2 = Inf;
chi2s = {};
labels = {};
fevals = [];
arWaitbar(0);
for j=1:length(filenames)
    arWaitbar(j,length(filenames));
    fname = ['./Results/' filenames{j} '/workspace.mat'];
    if(exist(fname,'file'))
        tmpple = load(fname);
        if(isfield(tmpple.ar, 'chi2s'))
            chi2s{end+1} = tmpple.ar.chi2s; %#ok<AGROW>
            labels{end+1} = filenames{j}; %#ok<AGROW>
            
            fevals(1:length(tmpple.ar.fun_evals),end+1) = tmpple.ar.fun_evals; %#ok<AGROW>
            
            minchi2 = min([minchi2 min(tmpple.ar.chi2s)]);
        end
    else
        fprintf('%s does not contains PLE\n', filenames{j});
    end
end
arWaitbar(-1);


figure(1)
h = nan(1,length(chi2s));
colors = jet(length(chi2s));
for j=1:length(chi2s)
    h(j) = semilogy(sort(chi2s{j}) + 1 - minchi2, 'o-', 'Color', colors(j,:), ...
        'MarkerFaceColor','w', ...
        'LineWidth',1, 'MarkerSize',4);
    hold on
end
hold off
legend(h, strrep(labels, '_', '\_'));
xlabel('run index (sorted by likelihood)');
ylabel('likelihood');

figure(2)
boxplot(log10(fevals), 'orientation', 'horizontal', 'labels', labels, ...
    'orientation', 'horizontal', ...
    'plotstyle', 'compact', 'colors', colors);
xlabel('log_{10} number of function evaluations');

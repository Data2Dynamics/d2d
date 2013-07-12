function arCompareFits(filenames, sortindex)

if(nargin==0)
    filenames = fileChooserMulti('./Results', true);
end
if(~exist('sortindex','var'))
    sortindex = -1;
end
if(~iscell(filenames))
    filelist = fileList('./Results');
    filenames = filelist(filenames);
end

minchi2 = Inf;
chi2s = {};
chi2sconstr = {};
optim_krit = {};
labels = {};
fevals = [];
timing = [];
exitflag = {};
arWaitbar(0);
jcount = 0;
for j=1:length(filenames)
    arWaitbar(j,length(filenames));
    fname = ['./Results/' filenames{j} '/workspace.mat'];
    if(exist(fname,'file'))
        tmpple = load(fname);
        if(isfield(tmpple.ar, 'chi2s'))
            jcount = jcount + 1;
            chi2s{jcount} = tmpple.ar.chi2s; %#ok<AGROW>
            chi2sconstr{jcount} = tmpple.ar.chi2sconstr; %#ok<AGROW>
            labels{jcount} = filenames{j}; %#ok<AGROW>
            if(isfield(tmpple.ar, 'optim_crit'))
                optim_krit{jcount} = tmpple.ar.optim_crit; %#ok<AGROW>
            end
            fevals(1:length(tmpple.ar.fun_evals),jcount) = tmpple.ar.fun_evals; %#ok<AGROW>
            timing(1:length(tmpple.ar.timing),jcount) = tmpple.ar.timing; %#ok<AGROW>
            exitflag{jcount} = tmpple.ar.exitflag; %#ok<AGROW>
            
            minchi2 = min([minchi2 min(tmpple.ar.chi2s)]);
        end
    else
        fprintf('%s does not contains PLE\n', filenames{j});
    end
end
arWaitbar(-1);

if(sortindex~=-1)
    [~, isort] = sort(chi2s{sortindex});
end

figure(1)
h = nan(1,length(chi2s));
colors = [];
for j=1:length(chi2s)
    if(sortindex==-1)
        [chi2s_sorted,isort] = sort(chi2s{j});
        exit_sorted = exitflag{j}(isort);
    else
        chi2s_sorted = chi2s{j}(isort);
        exit_sorted = exitflag{j}(isort);
    end
    if(~isempty(optim_krit{j}))
        optim_krit{j} = optim_krit{j}(isort); %#ok<AGROW>
    end
    C = arLineMarkersAndColors(j,length(chi2s),[],'no',[]);
    colors(j,:) = C{6}; %#ok<AGROW>
    h(j) = semilogy(chi2s_sorted + 1 - minchi2, '-', C{:});
    hold on
    C = arLineMarkersAndColors(j,length(chi2s),[],'o','none');
    semilogy(find(exit_sorted>0), chi2s_sorted(exit_sorted>0) + 1 - minchi2, C{:}, ...
        'MarkerFaceColor','w', ...
        'MarkerSize',4);
    C = arLineMarkersAndColors(j,length(chi2s),[],'x','none');
    semilogy(find(exit_sorted==0), chi2s_sorted(exit_sorted==0) + 1 - minchi2, C{:}, ...
        'MarkerSize',6);
end
hold off
legend(h, strrep(labels, '_', '\_'), 'Location','NorthWest');
xlabel('run index (sorted by likelihood)');
ylabel('likelihood');

figure(2)
subplot(3,1,1);
boxplot(log10(fevals), 'orientation', 'horizontal', 'labels', labels, ...
    'orientation', 'horizontal', ...
    'plotstyle', 'compact', 'colors', colors);
xlabel('log_{10} number of function evaluations');

subplot(3,1,2);
boxplot(log10(timing), 'orientation', 'horizontal', 'labels', labels, ...
    'orientation', 'horizontal', ...
    'plotstyle', 'compact', 'colors', colors);
xlabel('log_{10} runtime [s]');

subplot(3,1,3);
for j=1:length(chi2s)
    h(j) = semilogy(optim_krit{j}, 'o-', 'Color', colors(j,:), ...
        'MarkerFaceColor','w', ...
        'LineWidth',1, 'MarkerSize',4);
    hold on
end
hold off
% legend(h, strrep(labels, '_', '\_'));
xlabel('run index (sorted by likelihood)');
ylabel('first order optimality criterion');

if(~isempty(chi2sconstr))
    figure(3)
    h = nan(1,length(chi2sconstr));
    colors = [];
    for j=1:length(chi2sconstr)
        if(sortindex==-1)
            [~,isort] = sort(chi2s{j});
            chi2sconstr_sorted = chi2sconstr{j}(isort);
            exit_sorted = exitflag{j}(isort);
        else
            chi2sconstr_sorted = chi2sconstr{j}(isort);
            exit_sorted = exitflag{j}(isort);
        end
        C = arLineMarkersAndColors(j,length(chi2sconstr),[],'no',[]);
        colors(j,:) = C{6}; %#ok<AGROW>
        h(j) = semilogy(chi2sconstr_sorted, '-', C{:});
        hold on
        C = arLineMarkersAndColors(j,length(chi2s),[],'o','none');
        semilogy(find(exit_sorted>0), chi2sconstr_sorted(exit_sorted>0), C{:}, ...
            'MarkerFaceColor','w', ...
            'MarkerSize',4);
        C = arLineMarkersAndColors(j,length(chi2s),[],'x','none');
        semilogy(find(exit_sorted==0), chi2sconstr_sorted(exit_sorted==0), C{:}, ...
            'MarkerSize',6);
    end
    hold off
    legend(h, strrep(labels, '_', '\_'), 'Location','NorthWest');
    xlabel('run index (sorted by likelihood)');
    ylabel('constraint violation');
end
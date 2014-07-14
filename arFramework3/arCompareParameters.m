function arCompareParameters(filenames)

if(nargin==0)
    filenames = fileChooserMulti('./Results', true);
end
if(~iscell(filenames))
    filelist = fileList('./Results');
    filenames = filelist(filenames);
end

pLabel = {};
p = {};
pLabelCollect = {};
for j=1:length(filenames)
    fname = ['./Results/' filenames{j} '/workspace_pars_only.mat'];
    if(exist(fname,'file'))
        tmp = load(fname);
        pLabel{j} = tmp.ar.pLabel; %#ok<AGROW>
        pLabelCollect = union(pLabelCollect, pLabel{j});
        p{j} = tmp.ar.p; %#ok<AGROW>
    end
end

ps = nan(length(pLabelCollect),length(filenames));
for j=1:length(filenames)
    ps(ismember(pLabelCollect, pLabel{j}),j) = p{j};
end

figure(1)
for j=1:length(filenames)
    C = arLineMarkersAndColors(j, length(filenames), [], [], '-');
    h(j) = plot(ps(:,j), 1:length(pLabelCollect), C{:}); %#ok<AGROW>
    hold on
end
arSpacedAxisLimits;
for j=1:length(pLabelCollect)
    plot3(xlim, [j j], [-1 -1], '-', 'Color', [1 1 1]*0.8);
    hold on
end
plot3([0 0], ylim, [-0.5 -0.5], 'k-');
hold off

xlabel('parameter value');
set(gca, 'YTick', 1:length(pLabelCollect));
set(gca, 'YTickLabel', pLabelCollect);
ylim([0 length(pLabelCollect)+1])
legend(h, strrep(filenames,'_','\_'))

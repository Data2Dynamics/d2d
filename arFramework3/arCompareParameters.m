% arCompareParameters(filenames, onlyCommons)

function arCompareParameters(filenames, onlyCommons)

if(nargin==0 || isempty(filenames))
    filenames = fileChooserMulti('./Results', true);
end
if(~exist('onlyCommons','var'))
    onlyCommons = false;
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
        if(onlyCommons)
            if(~isempty(pLabelCollect))
                pLabelCollect = intersect(pLabelCollect, pLabel{j});
            else
                pLabelCollect = pLabel{j};
            end
        else
            pLabelCollect = union(pLabelCollect, pLabel{j});
        end
        p{j} = tmp.ar.p; %#ok<AGROW>
    end
end

ps = nan(length(pLabelCollect),length(filenames));
for j=1:length(filenames)
    if(onlyCommons)
        ps(:,j) = p{j}(ismember(pLabel{j}, pLabelCollect));
    else
        ps(ismember(pLabelCollect, pLabel{j}),j) = p{j};
    end
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
set(gca, 'YTickLabel', arNameTrafo(pLabelCollect));
ylim([0 length(pLabelCollect)+1])
legend(h, strrep(filenames,'_','\_'))

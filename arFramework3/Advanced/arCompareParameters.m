% arCompareParameters(filenames, onlyCommons, onlyFitted)
% 
% params = arCompareParameters
% 
% params = 
%            ps: [136x7 double]         parameters as in ar.p
%     filenames: {1x7 cell}             workspace names
%        pLabel: {136x1 cell}           parameter names
%        qLog10: [136x7 double]         log-scale? ar.qLog10
%        ps_log: [136x7 double]         all parameters at the log-scale


function params = arCompareParameters(filenames, onlyCommons, onlyFitted)

if(nargin==0 || isempty(filenames))
    filenames = fileChooserMulti('./Results', true);
end
if(~exist('onlyCommons','var'))
    onlyCommons = false;
end
if(~exist('onlyFitted','var'))
    onlyFitted = false;
end
if(~iscell(filenames))
    filelist = fileList('./Results');
    filenames = filelist(filenames);
end

pLabel = {};
p = {};
qlog = {};
pLabelCollect = {};
for j=1:length(filenames)
    fname = ['./Results/' filenames{j} '/workspace_pars_only.mat'];
    if(exist(fname,'file'))
        tmp = load(fname);
        if(onlyFitted)
            pLabel{j} = tmp.ar.pLabel(tmp.ar.qFit==1); %#ok<AGROW>
            p{j} = tmp.ar.p(tmp.ar.qFit==1); %#ok<AGROW>
            qlog{j} = tmp.ar.qLog10(tmp.ar.qFit==1); %#ok<AGROW>
        else
            pLabel{j} = tmp.ar.pLabel; %#ok<AGROW>
            p{j} = tmp.ar.p; %#ok<AGROW>
            qlog{j} = tmp.ar.qLog10; %#ok<AGROW>
        end
        if(onlyCommons)
            if(~isempty(pLabelCollect))
                pLabelCollect = intersect(pLabelCollect, pLabel{j});
            else
                pLabelCollect = pLabel{j};
            end
        else
            pLabelCollect = union(pLabelCollect, pLabel{j});
        end
    end
end

ps = nan(length(pLabelCollect),length(filenames));
qlogs = nan(length(pLabelCollect),length(filenames));
for j=1:length(filenames)
    if(onlyCommons)
        ps(:,j) = p{j}(ismember(pLabel{j}, pLabelCollect));
        qlogs(:,j) = qlog{j}(ismember(pLabel{j}, pLabelCollect));
    else
        ps(ismember(pLabelCollect, pLabel{j}),j) = p{j};
        qlogs(ismember(pLabelCollect, pLabel{j}),j) = qlog{j};
    end
end

if nargout>0
    params.ps = ps;
    params.filenames = filenames;
    params.pLabel = pLabelCollect;
    params.qLog10 = qlogs;
    params.ps_log = params.ps; % all parameters on the log-scale
    params.ps_log(params.qLog10==0) = log10(params.ps(params.qLog10==0));
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

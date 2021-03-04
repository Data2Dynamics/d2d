% arPlotPEtab(filename)
%
% Uses the information of a visualization*.tsv to plot data in ar.model.data.yExp
% to compare plots of PEtab models in different tools
%
% Example:
% arPlotPEtab('visualization_Crauste_CellSystems2017.tsv')
%
% Wiki: https://github.com/PEtab-dev/PEtab/blob/master/doc/documentation_data_format.rst#visualization-table

function arPlotPEtab(filename)

global ar;

if ~contains(filename,'.tsv')
    if ~contains(filename,'.')
        filename = [filename '.tsv'];
    else
        error('this file type is not supported!')
    end
end

%% Read in tsv file
T = tdfread(filename); % all data of the model
fn = fieldnames(T);
for i = 1:length(fn)
    if ischar(T.(fn{i}))
        T.(fn{i}) = regexprep(string(T.(fn{i})),' ','');
    end
end
%% Get plotIds
if isfield(T,'plotId')
    [~,~,idxplots] = unique(T.plotId);
else
    idxplots = 1:size(T.(fn{1}),1);
    warning('arLoadVisPETab.m: No plotId given. All plots are plotted in separate figures.')
end

%% Get infos for all rows with same plotId
for p=1:length(unique(idxplots))
    idxp = find(idxplots==p); 
    %% y
    if isfield(T,'yValues')
        [Obs,ia] = intersect(ar.model.data.y,T.yValues(idxp));
        y = ar.model.data.yExp(:,ia);
        if length(Obs)<length(T.yValues(idxp))
            warning(['arLoadVisPETab.m: Observable ' setdiff(T.yValues(idxp),Obs) ' not found in ar struct. Not plotted. Check name spelling.'])
        end
    end
    %% ystd
    if isfield(T,'plotTypeData')
        warning('arPlotPEtab: plotTypeData not implemented. In d2d each replicate is in the data. If noise is provided they are added as errorbars.')
    end
    if isfield(T,'datasetId')
        warning('arPlotPEtab: datasetId not yet implemented.')
    end
    if isfield(ar.model.data,'yExpStd') %&& ~all(isnan(ar.model.data.yExpStd))
        yStd = ar.model.data.yExpStd(:,ia);
    end

    %% x
    if ~isfield(T,'xValues') || (isfield(T,'xValues') && all(strcmp(T.xValues(idxp),'time')))
        if size(ar.model.data.tExp,2)==1
            x = repmat(ar.model.data.tExp,1,size(y,2));
        else
            x = ar.model.data.tExp(:,ia);
        end
    else
        % to do: dose curves
    end
    if isfield(T,'yOffset')
        y = y+T.yOffset(idxp)';
    end
    if isfield(T,'xOffset')
        x = x+T.xOffset(idxp)';
    end
    
    %% Scale
    if isfield(T,'xScale')
        if strcmp(T.xScale,'log')
            x = log(x);
        elseif strcmp(T.xScale,'log10')
            x = log10(x);
        end
    end
    if isfield(T,'yScale')
        if strcmp(T.yScale,'log')
            y = log(y);
        elseif strcmp(T.yScale,'log10')
            y = log10(y);
        end
    end
    
    %% Figure
    figure;
    if isfield(T,'plotTypeSimulation') && all(strcmp(T.plotTypeSimulation(idxp),'ScatterPlot'))
        hold on
        for i=1:size(y,2)
            scatter(x(:,i),y(:,i),'filled')
        end
        if exist('yStd','var')
            errorbar(x,y,yStd,'LineStyle','none');
        end        
    elseif isfield(T,'plotTypeSimulation') && all(strcmp(T.plotTypeSimulation(idxp),'BarPlot'))
        bar(x,y)
        if exist('yStd','var')
            hold on
            errorbar(x,y,y-yStd,y+yStd,'Color',[0 0 0],'LineStyle','none');
        end
    else
        if exist('yStd','var')
            errorbar(x,y,yStd,'.-','MarkerSize',10)
        else
            plot(x,y,'.-','MarkerSize',10)
        end
    end
    if isfield(T,'xLabel')
        xlabel(T.xLabel(idxp(1)))
    end
    if isfield(T,'yLabel')
        ylabel(T.yLabel(idxp(1)))
    end
    if isfield(T,'legendEntry')
        legend(T.legendEntry(idxp))
    end
    if isfield(T,'plotName')
        title(T.plotName(idxp(1)))
    end
end
    


%ar.model.qPlotYs([1, 2, 8, 14]) = 1; 

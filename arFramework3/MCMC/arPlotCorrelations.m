% function arPlotCorrelations(ParaPerFigure,ParaVector)
%
%  NumberOfPlottedParametersPerFigure
%            Number of parameters that are plotted into same figure 
%
%  ParaVector
%       Give parameter vector to only plot correlations of certain
%       parameters, e.g. ParaVector = [ 1 2 3 4 6 9 112 114]
%
%

function arPlotCorrelations(varargin)

if nargin > 1
    ParaVector = varargin{2};    
else 
    ParaVector = find(ar.qFit==1);
end

if nargin > 0
    ParaPerFigure = varargin{1};    
else
    ParaPerFigure = 10;
end

global ar


NumberOfFigures = ceil(length(ParaVector)/ParaPerFigure);

for FigureNumber = 1:NumberOfFigures

    if FigureNumber < NumberOfFigures    
        FigureParaVector = ParaVector(1+(FigureNumber-1)*ParaPerFigure:(FigureNumber-1)*ParaPerFigure+ParaPerFigure);
    else
        FigureParaVector = ParaVector(1+(FigureNumber-1)*ParaPerFigure:end);
    end
    NumberOfColumns = length(FigureParaVector)+1;
    NumberOfRows = length(FigureParaVector);

    figure('Name',sprintf('2D Correlations - Figure %d / %d',FigureNumber,NumberOfFigures),'NumberTitle','off', 'Position', [10 300 2000 1000])
    % Plot parameter distribution
    for paraindex=FigureParaVector
        hold on
        i = find(paraindex==FigureParaVector);
        if i ~=1
            subplot(NumberOfRows,NumberOfColumns,1+(i-1)*NumberOfColumns)
            histogram(ar.ps(:,paraindex))
            line([ar.p(paraindex), ar.p(paraindex)], ylim, 'LineWidth', 2, 'Color', 'r');
            title(ar.pLabel{paraindex});
            ax{i} = gca; % current axes
            camroll(90)
        end
        subplot(NumberOfRows,NumberOfColumns,1+i+(i-1)*NumberOfColumns)
        histogram(ar.ps(:,paraindex))
        line([ar.p(paraindex), ar.p(paraindex)], ylim, 'LineWidth', 2, 'Color', 'r');
        title(ar.pLabel{paraindex});
        ax{i} = gca; % current axes
    % Plot 2D-correlations
        for CorrelationIndex=FigureParaVector(1:i-1)
            k = find(CorrelationIndex==FigureParaVector);
            R = corrcoef(ar.ps(:,paraindex,1),ar.ps(:,CorrelationIndex,1));
            hold on
            subplot(NumberOfRows,NumberOfColumns,1+k+(i-1)*NumberOfColumns)
            [N,C]= hist3([ar.ps(:,paraindex),ar.ps(:,CorrelationIndex)],[80 80]);
            imagesc([C{2}(1),C{2}(80)],[C{1}(1),C{1}(80)],log(N));
            set(gca,'YDir','normal')
            hold on             
            title(sprintf('R = %0.2f',R(2)));


           xlim(ax{k}.XLim);
           ylim(ax{i}.XLim);

        end
    end      
end

end

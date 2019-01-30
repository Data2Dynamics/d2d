% [N, centers] = arPlotHist(data, edges, xf, yf, vertical)
% Plots Histogram and shows scaling of empirical data xf,yf


function [N, centers] = arPlotHist(data, edges, xf, yf, vertical)

if(~exist('vertical','var'))
    vertical = false;
end

N = histc(data,edges);
N = [N(1:end-2); N(end-1)+N(end)];

centers = (edges(1:end-1) + edges(2:end))/2;
if(~vertical)
    bar(centers,N,1,'FaceColor',[.8 .8 .8]);
else
    barh(centers,N,1,'FaceColor',[.8 .8 .8]);
end

if(nargin>2 && ~isempty(xf) && ~isempty(yf))
    qh = ishold;
    if(~qh)
        hold on
    end
    scale = arMatchHistWithEmpFunction(centers,N,xf,yf);
    if(~vertical)
        plot(xf, scale*yf, 'r-', 'LineWidth',1);
    else
        plot(scale*yf, xf, 'r-', 'LineWidth',1);
    end
    if(qh)
        hold on
    else
        hold off
    end
end

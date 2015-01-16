% plot sampling of likelihood

function arPlotSample

global ar
global pleGlobals

figure(1)
clf;

% constants
overplot = 0.1;

if(ar.config.fiterrors == 1)
    chi2s = 2*ar.ndata*log(sqrt(2*pi)) + ar.sampling.chi2s;
    chi2curr = 2*ar.ndata*log(sqrt(2*pi)) + ar.chi2fit;
    ylabeltmp = '-2*log(L)';
else
    chi2s = ar.sampling.chi2s;
    chi2curr = ar.chi2fit;
    ylabeltmp = '\chi^2';
end

if(length(ar.sampling.ps)==1)
    plot(ar.sampling.ps{1}, chi2s, 'k');
    hold on
    plot(ar.p(ar.sampling.index), chi2curr, '*k');
    if(isfield(ar.sampling, 'chi2s_ple') && ...
            length(ar.sampling.ps{1})==length(ar.sampling.chi2s_ple))
        plot(ar.sampling.ps{1}, ar.sampling.chi2s_ple, 'r--');
        qerrors = ar.sampling.ple_errors == 0;
        plot(ar.sampling.ps{1}(qerrors), ar.sampling.chi2s_ple(qerrors), 'ro');
        qerrors = ar.sampling.ple_errors < 0;
        plot(ar.sampling.ps{1}(qerrors), ar.sampling.chi2s_ple(qerrors), 'r*');
        qerrors = ar.sampling.ple_errors > 1;
        plot(ar.sampling.ps{1}(qerrors), ar.sampling.chi2s_ple(qerrors), 'rs');
    end
    hold off
    spacedAxisLimits(gca, overplot)
    ylabel(ylabeltmp);
    xlabel(arNameTrafo(ar.pLabel{ar.sampling.index}));
elseif(length(ar.sampling.ps)==2)
    [X,Y] = meshgrid(ar.sampling.ps{1}, ar.sampling.ps{2});
    
    dchi2 = chi2inv(0.95, 1);
    Cmin = min(chi2s(:));
    
    Cmax = Cmin+dchi2*5;
    colormap(gray); 
%     colormap(jet);
    
%     Cvals = logspace(log10(Cmin), log10(Cmax), 100);
%     contourf(X,Y,chi2s, Cvals, 'LineColor', 'none');

%     Cvals = logspace(log10(Cmin), log10(Cmax), 20);
    Cvals = linspace(Cmin, Cmax, 30);
    contour(X,Y,chi2s, Cvals);
    
    h = colorbar('EastOutside');
    title(h, ylabeltmp);
    hold(h, 'on');
    plot(h, [-5 5], [0 0]+Cmin+dchi2, 'r', 'LineWidth', 2);
    hold(h, 'off');
    
    hold on
    if(isfield(ar, 'ps') && ~isempty(ar.ps))
        plot(ar.ps(:,ar.sampling.index(1)), ar.ps(:,ar.sampling.index(2)), 'xb', 'MarkerSize', 4);
    end
    contour(X,Y,chi2s, Cmin+dchi2,'r', 'LineWidth', 1)
%     plot(ar.p(ar.sampling.index(1)), ar.p(ar.sampling.index(2)), '*w');
    plot(ar.p(ar.sampling.index(1)), ar.p(ar.sampling.index(2)), '*k');
    
    if(isfield(ar,'covar'))
        error_ellipse(ar.covar([ar.sampling.index(1) ar.sampling.index(2)], ...
            [ar.sampling.index(1) ar.sampling.index(2)]), ...
            [ar.p(ar.sampling.index(1)) ar.p(ar.sampling.index(2))], 'style', 'k--', 'conf', 0.95)
    end
    
    if(~isempty(pleGlobals) && ~isempty(pleGlobals.ps{ar.sampling.index(1)}))
        plot(pleGlobals.ps{ar.sampling.index(1)}(:,ar.sampling.index(1)), ...
            pleGlobals.ps{ar.sampling.index(1)}(:,ar.sampling.index(2)), ...
            'k--', 'LineWidth', 1);
    end
    if(~isempty(pleGlobals) && ~isempty(pleGlobals.ps{ar.sampling.index(2)}))
        plot(pleGlobals.ps{ar.sampling.index(2)}(:,ar.sampling.index(1)), ...
            pleGlobals.ps{ar.sampling.index(2)}(:,ar.sampling.index(2)), ...
            'k--', 'LineWidth', 1);
    end
    
    hold off
    
    xlabel(arNameTrafo(ar.pLabel{ar.sampling.index(1)}));
    ylabel(arNameTrafo(ar.pLabel{ar.sampling.index(2)}));
    
    xlim([min(ar.sampling.ps{1}) max(ar.sampling.ps{1})])
    ylim([min(ar.sampling.ps{2}) max(ar.sampling.ps{2})])
    
end


function spacedAxisLimits(g, overplot)
[xmin xmax ymin ymax] = axisLimits(g);
xrange = xmax - xmin;
if(xrange == 0)
	xrange = 1;
end
yrange = ymax - ymin;
if(yrange == 0)
	yrange = 1;
end
xlim(g, [xmin-(xrange*overplot) xmax+(xrange*overplot)]);
ylim(g, [ymin-(yrange*overplot) ymax+(yrange*overplot)]);



function [xmin xmax ymin ymax] = axisLimits(g)
p = get(g,'Children');
xmin = nan;
xmax = nan;
ymin = nan;
ymax = nan;
for j = 1:length(p)
	if(~strcmp(get(p(j), 'Type'), 'text'))
		xmin = min([xmin toRowVector(get(p(j), 'XData'))]);
		xmax = max([xmax toRowVector(get(p(j), 'XData'))]);
		%         get(p(j), 'UData')
		%         get(p(j), 'LData')
		%         set(p(j), 'LData', get(p(j), 'LData')*2)
		if(strcmp(get(p(j), 'Type'),'hggroup'))
			ymin = min([ymin toRowVector(get(p(j), 'YData'))-toRowVector(get(p(j), 'LData'))]);
			ymax = max([ymax toRowVector(get(p(j), 'YData'))+toRowVector(get(p(j), 'UData'))]);
		else
			ymin = min([ymin toRowVector(get(p(j), 'YData'))]);
			ymax = max([ymax toRowVector(get(p(j), 'YData'))]);
		end
	end
end



function b = toRowVector(a)
b = a(:)';


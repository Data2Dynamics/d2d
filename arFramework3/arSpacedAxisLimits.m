function arSpacedAxisLimits(g, overplot)
if(~exist('g','var'))
    g = gca;
end
if(~exist('overplot','var'))
    overplot = 0.1;
end
[xmin xmax ymin ymax] = axisLimits(g);

% x-axis
if(~strcmp(get(g, 'XScale'), 'linear'))
    xmin = log10(xmin);
    xmax = log10(xmax);
end
xrange = xmax - xmin;
if(xrange == 0)
    xrange = 1;
end
if(~isnan(xrange))
    if(strcmp(get(g, 'XScale'), 'linear'))
        xlim(g, [xmin-(xrange*overplot) xmax+(xrange*overplot)]);
    else
        xlim(g, 10.^[xmin-(xrange*overplot) xmax+(xrange*overplot)]);
    end
end

% y-axis
if(~strcmp(get(g, 'YScale'), 'linear'))
    ymin = log10(ymin);
    ymax = log10(ymax);
end
yrange = ymax - ymin;
if(yrange == 0)
    yrange = 1;
end
if(~isnan(yrange))
    if(strcmp(get(g, 'YScale'), 'linear'))
        ylim(g, [ymin-(yrange*overplot) ymax+(yrange*overplot)]);
    else
        ylim(g, 10.^[ymin-(yrange*overplot) ymax+(yrange*overplot)]);
    end
end


function [xmin xmax ymin ymax] = axisLimits(g)
p = get(g,'Children');
xmin = nan;
xmax = nan;
ymin = nan;
ymax = nan;
for j = 1:length(p)
	if(~strcmp(get(p(j), 'Type'), 'text'))
        xtmp = toRowVector(get(p(j), 'XData'));
        xtmp = xtmp(~isinf(xtmp));
		xmin = min([xmin xtmp]);
		xmax = max([xmax xtmp]);
        
        ytmp = toRowVector(get(p(j), 'YData'));
        ytmp = ytmp(~isinf(ytmp));
		if(strcmp(get(p(j), 'Type'),'hggroup'))
            ytmpu = toRowVector(get(p(j), 'LData'));
            ytmpl = toRowVector(get(p(j), 'UData'));
            ytmpu = ytmpu(~isinf(ytmpu));
            ytmpl = ytmpl(~isinf(ytmpl));
        
			ymin = min([ymin ytmp-ytmpu]);
			ymax = max([ymax ytmp+ytmpl]);
		else
			ymin = min([ymin ytmp]);
			ymax = max([ymax ytmp]);
		end
	end
end

function b = toRowVector(a)
b = a(:)';



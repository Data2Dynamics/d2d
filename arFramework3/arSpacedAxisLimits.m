function arSpacedAxisLimits(g, overplot)
if(~exist('g','var'))
    g = gca;
end
if(~exist('overplot','var'))
    overplot = 0.1;
end
[xmin xmax ymin ymax] = axisLimits(g);
xrange = xmax - xmin;
if(xrange == 0)
    xrange = 1;
end
yrange = ymax - ymin;
if(yrange == 0)
    yrange = 1;
end
if(~isnan(xrange))
    xlim(g, [xmin-(xrange*overplot) xmax+(xrange*overplot)]);
end
if(~isnan(yrange))
    ylim(g, [ymin-(yrange*overplot) ymax+(yrange*overplot)]);
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



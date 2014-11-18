function arSpacedAxisLimits(g, overplot)
if(~exist('g','var'))
    g = gca;
end
if(~exist('overplot','var'))
    overplot = 0.1;
end

if(length(g)==1)
    [xmin, xmax, ymin, ymax] = axisLimits(g);
else
    xmin = Inf;
    xmax = -Inf;
    ymin = Inf;
    ymax = -Inf;
    for j=1:length(g)
        [tmpxmin, tmpxmax, tmpymin, tmpymax] = axisLimits(g(j));
        if(tmpxmin<xmin)
            xmin = tmpxmin;
        end
        if(tmpxmax>xmax)
            xmax = tmpxmax;
        end
        if(tmpymin<ymin)
            ymin = tmpymin;
        end
        if(tmpymax>ymax)
            ymax = tmpymax;
        end
    end
end

for j=1:length(g)
    % x-axis
    if(~strcmp(get(g(j), 'XScale'), 'linear'))
        xmin = log10(xmin);
        xmax = log10(xmax);
    end
    xrange = xmax - xmin;
    if(xrange == 0)
        xrange = 1;
    end
    if(~isnan(xrange))
        if(strcmp(get(g(j), 'XScale'), 'linear'))
            xlim(g(j), [xmin-(xrange*overplot) xmax+(xrange*overplot)]);
        else
            xlim(g(j), 10.^[xmin-(xrange*overplot) xmax+(xrange*overplot)]);
        end
    end
    
    % y-axis
    if(~strcmp(get(g(j), 'YScale'), 'linear'))
        ymin = log10(ymin);
        ymax = log10(ymax);
    end
    yrange = ymax - ymin;
    if(yrange == 0)
        yrange = 1;
%     if(yrange < 1e-3)
%         yrange = 1e-3;
    end
    if(~isnan(yrange))
        if(strcmp(get(g(j), 'YScale'), 'linear'))
            if(ymin-(yrange*overplot) < ymax+(yrange*overplot))
                ylim(g(j), [ymin-(yrange*overplot) ymax+(yrange*overplot)]);
            end
        else
            if(ymin-(yrange*overplot) < ymax+(yrange*overplot))
                ylim(g(j), 10.^[ymin-(yrange*overplot) ymax+(yrange*overplot)]);
            end
        end
    end
end

function [xmin, xmax, ymin, ymax] = axisLimits(g)
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
        end
        ymin = min([ymin ytmp]);
        ymax = max([ymax ytmp]);
    end
end

function b = toRowVector(a)
b = a(:)';



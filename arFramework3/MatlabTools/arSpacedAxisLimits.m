% arSpacedAxisLimits(g, overplot, doX, doY, minYrange)
%
% g: empty -> gca
%    axis handle
%    vector of axis handles
%    figure handle -> get children axis handles
%
% overplot: how much space around elements in figure in % (default is 0.1)
%
% doX: boolean -> consider x-axis? (default is true)
%      vector -> vector for [xmin xmax], overrides overplot functionality
%
% doY: boolean -> consider y-axis? (default is true)
%      vector -> vector for [ymin ymax], overrides overplot functionality
% 
% minYrange: minimum y-range, overrides overplot functionality

function arSpacedAxisLimits(g, overplot, doX, doY, minYrange)
if(~exist('g','var'))
    g = gca;
end
if(isa(g,'matlab.ui.Figure'))
    h = g;
    g = [];
    for j=1:length(h.Children)
        if(isa(h.Children(j),'matlab.graphics.axis.Axes'))
            g(end+1) = h.Children(j); %#ok<AGROW>
        end
    end
end
if(~exist('overplot','var') || isempty(overplot))
    overplot = 0.1;
end
if(~exist('doX','var') || isempty(doX))
    doX = true;
end
if(~exist('doY','var') || isempty(doY))
    doY = true;
end

overplotx = overplot;
overploty = overplot;

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

if(length(doX)>1)
    xmin = doX(1);
    xmax = doX(2);
    doX = true;
    overplotx = 0;
end
if(length(doY)>1)
    ymin = doY(1);
    ymax = doY(2);
    doY = true;
    overploty = 0;
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
    if(~isnan(xrange) && doX)
        if(strcmp(get(g(j), 'XScale'), 'linear'))
            xlim(g(j), [xmin-(xrange*overplotx) xmax+(xrange*overplotx)]);
        else
            xlim(g(j), 10.^[xmin-(xrange*overplotx) xmax+(xrange*overplotx)]);
        end
    end
    
    % y-axis
    if(~strcmp(get(g(j), 'YScale'), 'linear'))
        ymin = log10(ymin);
        ymax = log10(ymax);
    end
    yrange = ymax - ymin;
    if(exist('minYrange','var') && yrange < minYrange)
        yrange = minYrange;
    end
    if(yrange == 0)
        yrange = 1;
    end
    if(~isnan(yrange) && doY)
        if(strcmp(get(g(j), 'YScale'), 'linear'))
            if(ymin-(yrange*overploty) < ymax+(yrange*overploty))
                ylim(g(j), [ymin-(yrange*overploty) ymax+(yrange*overploty)]);
            end
        else
            if(ymin-(yrange*overploty) < ymax+(yrange*overploty))
                ylim(g(j), 10.^[ymin-(yrange*overploty) ymax+(yrange*overploty)]);
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



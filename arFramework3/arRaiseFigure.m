function [h, fastPlotTmp] = arRaiseFigure(base, fieldname, figname, ...
    figcount, fastPlot, level)

if(~exist('level','var'))
    level = 1;
end
if(~exist('fastPlot','var'))
    fastPlot = false;
end
if(~exist('figcount','var'))
    figcount = 1;
end

levels = [0.1 0.4 0.8];

openfigs = get(0,'Children');

figcolor = [1 1 1];
figdist = 0.02;

if(nargin<4)
    figpos = [0.1 0.1];
    figsize = [0.6 0.8];
else
    figpos = [levels(level)+((figcount-1)*figdist) 0.35-((figcount-1)*figdist)];
    figsize = [0.3 0.45];
end

base.time = now;
fastPlotTmp = fastPlot;

if(isfield(base, fieldname) && ~isempty(base.(fieldname)) && ...
        base.(fieldname) ~= 0 && ...
        sum(base.(fieldname)==openfigs)>0 && ...
        strcmp(get(base.(fieldname), 'Name'), figname))
    
    h = base.(fieldname);
    if(~fastPlot)
        figure(h);
    end
else
    h = figure('Name', figname, 'NumberTitle','off', ...
        'Units', 'normalized', 'Position', ...
        [figpos figsize]);
    set(h,'Color', figcolor);
    base.(fieldname) = h;
    fastPlotTmp = false;
end

if(~fastPlot)
    clf
end

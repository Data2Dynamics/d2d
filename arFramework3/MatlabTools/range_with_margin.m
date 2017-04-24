function [lb, ub] = range_with_margin(data, overplot)

if(~exist('overplot','var') || isempty(overplot))
    overplot = 0.1;
end

xmax = max(data);
xmin = min(data);

xrange = xmax - xmin;
if(xrange == 0)
    xrange = 1;
end

lb = xmin-(xrange*overplot);
ub = xmax+(xrange*overplot);
    
function [bounds,boundsNegative] = DefaultLbUbTransient(im,id)
global ar
if(nargin==0)
    im = 1;
    id = 1;
end

y = ar.model(im).data(id).yExp(~isnan(ar.model(im).data(id).yExp));
t = ar.model(im).data(id).tExp;
t = unique(t(~isnan(t)));
D = range(y);

bounds.pLabel = ar.pLabel;
bounds.lb = ar.lb;
bounds.ub = ar.ub;

ind = ar.fit_transient.indp.amp_sust;
if(~isempty(ind))
    bounds.lb(ind) = 1e-10;
    bounds.ub(ind) = D*2;
    if(ar.qLog10(ind)==1)
        bounds.lb(ind) = log10(bounds.lb(ind));
        bounds.ub(ind) = log10(bounds.ub(ind));
    end
end

ind = ar.fit_transient.indp.amp_trans;
if(~isempty(ind))
    bounds.lb(ind) = 1e-10;
    bounds.ub(ind) = D*2;
    if(ar.qLog10(ind)==1)
        bounds.lb(ind) = log10(bounds.lb(ind));
        bounds.ub(ind) = log10(bounds.ub(ind));
    end
end

indOffset = ar.fit_transient.indp.offset;
if(~isempty(indOffset))
    bounds.lb(indOffset) = nanmin(y)-0.5*D;
    bounds.ub(indOffset) = nanmax(y);
end

ind = ar.fit_transient.indp.timescale_sust;
if(~isempty(ind))
    bounds.lb(ind) = min(diff(t(:)))/2; %(t(2)-t(1))/2;
    bounds.ub(ind) = 2*(t(end)-t(1));
    if(ar.qLog10(ind)==1)
        bounds.lb(ind) = log10(bounds.lb(ind));
        bounds.ub(ind) = log10(bounds.ub(ind));
    end
end

ind = ar.fit_transient.indp.timescale_trans;
if(~isempty(ind))
    bounds.lb(ind) = (t(2)-t(1))/2;
    bounds.ub(ind) = 2*(t(end)-t(1));
    if(ar.qLog10(ind)==1)
        bounds.lb(ind) = log10(bounds.lb(ind));
        bounds.ub(ind) = log10(bounds.ub(ind));
    end
end


ind = ar.fit_transient.indp.sd;
if(~isempty(ind))
    bounds.ub(ind) = std(y)./sqrt(length(y)./length(t));
    if(ar.qLog10(ind)==1)
        bounds.ub(ind) = log10(bounds.ub(ind));
    end
end

%% for fits in negative direction:
boundsNegative = bounds;

if(~isempty(indOffset))
    boundsNegative.lb(indOffset) = nanmin(y);
    boundsNegative.ub(indOffset) = nanmax(y)+0.5*D;
end

    
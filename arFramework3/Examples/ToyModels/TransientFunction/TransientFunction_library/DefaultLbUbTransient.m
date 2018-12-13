% Pay attention: Defintion of bounds depends on the data.


function [bounds,boundsNegative] = DefaultLbUbTransient
global ar

bounds.pLabel = ar.pLabel;
bounds.lb = ar.lb;
bounds.ub = ar.ub;

%% The transient fuctions looks like this:
%     y1 = p4*(1-exp(-t/timescale_sust)).*exp(-p2*t);
%     y2 = p3*(1-exp(-t/timescale_sust));
%     y = y1+y2+p(5);
% 
%     bound.lb = [.5/(t(end)-t(1)),.5/(t(end)-t(1)),1e-10,1e-10,  nanmin(data)-0.5*D];
%     bound.ub = [2/(t(2)-t(1)),2/(t(2)-t(1)),D*2,D*2,  nanmax(data)];

%% amplitudes
ind = ar.fit_transient.indp.amp_sust;
for i=1:length(ind)
    [t,y,D] = getData(ar.pLabel{ind(i)});
    bounds.lb(ind(i)) = D*1e-6;
    bounds.ub(ind(i)) = D*2;
    if(ar.qLog10(ind(i))==1)
        bounds.lb(ind(i)) = log10(bounds.lb(ind(i)));
        bounds.ub(ind(i)) = log10(bounds.ub(ind(i)));
    end
end

ind = ar.fit_transient.indp.amp_trans;
for i=1:length(ind)
    [t,y,D] = getData(ar.pLabel{ind(i)});
    bounds.lb(ind(i)) = D*1e-6;
%     bounds.ub(ind(i)) = D*4; % this choice might be better for more dense data
    bounds.ub(ind(i)) = D*2; % this choice is more conservatvie and might be better for roughly sampled data.
    if(ar.qLog10(ind(i))==1)
        bounds.lb(ind(i)) = log10(bounds.lb(ind(i)));
        bounds.ub(ind(i)) = log10(bounds.ub(ind(i)));
    end
end

%% offsets
indOffset = ar.fit_transient.indp.offset;
for i=1:length(indOffset)
    [t,y,D] = getData(ar.pLabel{indOffset(i)});
    bounds.lb(indOffset(i)) = nanmin(y)-0.5*D;
    bounds.ub(indOffset(i)) = nanmax(y);
end
if(sum(ar.qLog10(indOffset(i)))>0)
    warning('Negative offsets are often required. Please set ar.qLog10 for offset parameters =0.');
end

%% timescales
ind = ar.fit_transient.indp.timescale_sust;
for i=1:length(ind)
    [t,y,D] = getData(ar.pLabel{ind(i)});
    if length(t)<10
        bounds.lb(ind(i)) = (t(2)-t(1))/2;
    else
        bounds.lb(ind(i)) = (t(2)-t(1));
    end
    bounds.ub(ind(i)) = 2*(t(end)-t(1));
    if(ar.qLog10(ind(i))==1)
        bounds.lb(ind(i)) = log10(bounds.lb(ind(i)));
        bounds.ub(ind(i)) = log10(bounds.ub(ind(i)));
    end
end

ind = ar.fit_transient.indp.timescale_trans;
for i=1:length(ind)
    [t,y,D] = getData(ar.pLabel{ind(i)});
    bounds.lb(ind(i)) = (t(2)-t(1))/2;
    bounds.ub(ind(i)) = 2*(t(end)-t(1));
    if(ar.qLog10(ind(i))==1)
        bounds.lb(ind(i)) = log10(bounds.lb(ind(i)));
        bounds.ub(ind(i)) = log10(bounds.ub(ind(i)));
    end
end

%% SDs
ind = ar.fit_transient.indp.sd;
for i=1:length(ind)
    [t,y,D] = getData(ar.pLabel{ind(i)});

    % lower bound:
    R = range(y);
    if R>0 && ~isnan(R) && ~isinf(R)
        bounds.lb(ind(i)) = R / 1e4;
        if(ar.qLog10(ind(i))==1)
            bounds.lb(ind(i)) = log10(bounds.lb(ind(i)));
        end
    else % keep bounds from ar.lb
    end
    
    % upper bound:
    SD = nanstd(y);
    % handle case SD=0 or SD too small:
    sdlb = bounds.lb(ind(i));
    if(ar.qLog10(ind(i))==1)
        sdlb = 10.^sdlb;
    end
    if isnan(SD) || SD<sdlb
        SD = sdlb*10;
    end
        
    bounds.ub(ind(i)) = SD./sqrt(length(y)./length(t));
    
    if(ar.qLog10(ind(i))==1)
        bounds.ub(ind(i)) = log10(bounds.ub(ind(i)));
    end
end

%% for fits in negative direction:
boundsNegative = bounds;

if(~isempty(indOffset))
    boundsNegative.lb(indOffset) = nanmin(y);
    boundsNegative.ub(indOffset) = nanmax(y)+0.5*D;
end

function [t,y,D] = getData(pname)
global ar
y = ar.fit_transient.(pname).yexp;
t = ar.fit_transient.(pname).t;

t = unique(t);
D = range(y);

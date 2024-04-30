% This function calls arSetPars but also translates info (e.g. bounds) into
% ar.fit_transient
function arSetParsTransient(pLabel, varargin)
global ar
if length(varargin)>0
    p = varargin{1};
else
    p = [];
end
if length(varargin)>1
    qFit = varargin{2};
else
    qFit = [];
end
if length(varargin)>2
    qLog10 = varargin{3};
else
    qLog10 = [];
end
if length(varargin)>3
    lb = varargin{4};
else
    lb = [];
end
if length(varargin)>4
    ub = varargin{5};
else
    ub = [];
end

% translate bounds to fit_transient.bounds
if(isfield(ar,'fit_transient'))
    ind = strmatch(pLabel,ar.fit_transient.bounds.pLabel,'exact');
    if ~isemtpy(lb)
        ar.fit_transient.bounds.lb(ind) = lb;
    end
    if ~isemtpy(ub)
        ar.fit_transient.bounds.ub(ind) = ub;
    end
end
% arSetPars(pLabel, [p], [qFit], [qLog10], [lb], [ub], [type], [meanp], [stdp])
arSetPars(pLabel_or_ar, varargin{:})
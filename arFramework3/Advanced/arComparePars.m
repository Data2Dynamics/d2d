
function pass = arComparePars(ar1,ar2,silent,bounds,fitted)

if exist('silent','var') && silent
    warning('off');
end

% Tolerances for a pass
rtol = 1e-3;
atol = 1e-4;

% Sort parameters of ar2 according to pLabels in ar1
p = nan(1,length(ar1.pLabel));
qlog = nan(1,length(ar1.pLabel));
for j=1:length(ar2.p)
    p(ismember(ar1.pLabel, ar2.pLabel{j})) = ar2.p(j);
    qlog(ismember(ar1.pLabel, ar2.pLabel{j})) = ar2.qLog10(j);
end

% if not present, give warning
if any(isnan(p))
    warning(['arComparePars.m: Parameter ' num2str(find(isnan(p))) ' (' ar2.pLabel{isnan(p)} ') are not present in SBML model.\n'])
end

% check correct log
checklog = ar1.qLog10(~isnan(qlog)) == qlog(~isnan(qlog));
if any(checklog~=1)
    if ar1.qLog10(checklog~=1)==1
        p(checklog~=1) = log10(p(checklog~=1));
        ar2.lb(checklog~=1) = log10(ar2.lb(checklog~=1));
        ar2.ub(checklog~=1) = log10(ar2.ub(checklog~=1));
    elseif qlog(checklog~=1)==1
        p(checklog~=1) = 10.^(p(checklog~=1));
        ar2.lb(checklog~=1) = 10.^(ar2.lb(checklog~=1));
        ar2.ub(checklog~=1) = 10.^(ar2.ub(checklog~=1));
    end
    warning('arComparePars.m: ar.qLog10 are different. SBML parameter and bounds are transformed for comparing.');
end

% Compare parameters
p1    = bsxfun(@max, ar1.p, atol);
p2    = bsxfun(@max, p, atol);
maxDifference = nanmax( ( (p1 - p2) ./ p2 ).^2 );

% Parameters acceptable?
if ( max( maxDifference ) > rtol )
    pass = 0;
else
    pass = 1;
end

% Same parameters fitted?
if ~exist('fitted','var') || fitted
    if any(ar1.qFit ~= ar2.qFit)
        warning(['arComparePars.m: qFit for parameters ' ar1.pLabel{ar1.qFit~=ar2.qFit} ' are different. Check it!']);
        pass = 0;
    end
end

% Same bounds?
if ~exist('bounds','var') || bounds
    if any(ar1.lb ~= ar2.lb)
        warning(['arComparePars.m: Lower bound for parameters ' num2str(find(ar1.lb~=ar2.lb)) ' (' ar1.pLabel{ar1.lb~=ar2.lb} ') are different. Check it!']);
        pass = 0;
    end 
    if any(ar1.ub ~= ar2.ub)
        warning(['arComparePars.m: Upper bound for parameters ' num2str(find(ar1.ub~=ar2.ub)) ' (' ar1.pLabel{ar1.ub~=ar2.ub} ') are different. Check it!']);
        pass = 0;
    end 
end
warning('on')


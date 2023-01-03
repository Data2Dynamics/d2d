
function arSetParsBounds(range)

global ar

if ~exist('range','var') || isempty(range)
    range = 2;
end
for i=1:length(ar.p)
    if ar.qLog10(i)
        ar.lb(i) = ar.p(i)-range;
        ar.ub(i) = ar.p(i)+range;
    elseif ar.p(i)>0
        ar.lb(i) = 10.^(log10(abs(ar.p(i)))-range);
        ar.ub(i) = 10.^(log10(abs(ar.p(i)))+range);
    elseif ar.p(i)<0
        ar.lb(i) = -10.^(log10(abs(ar.p(i)))+range);
        ar.ub(i) = -10.^(log10(abs(ar.p(i)))-range);
    else % if ar.p=0 und qLog10=0
        ar.lb(i) = -range;
        ar.ub(i) = range;
    end
end
function c = countElements(d, sorted)

if(~exist('sorted','var'))
    sorted = false;
end

if(istable(d))
    d = table2cell(d);
end

[ud, ~, jd] = unique(d);

maxlength = max(cellfun(@length, d));

counts = zeros(1,length(ud));
for j=1:length(ud)
    counts(j) = sum(jd==j);
end

if(sorted)
    [~, isort] = sort(counts);
    counts = counts(isort);
    ud = ud(isort);
end
if nargout==0
    for j=1:length(ud)
        fprintf([sprintf('%%%is', maxlength) ' : %i\n'], ud{j}, counts(j));
    end
else
    c = counts;
end

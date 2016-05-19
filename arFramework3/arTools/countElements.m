function countElements(d)

[ud, ~, jd] = unique(d);

for j=1:length(ud)
    fprintf('%s\t: #%i\n', ud{j}, sum(jd==j));
end
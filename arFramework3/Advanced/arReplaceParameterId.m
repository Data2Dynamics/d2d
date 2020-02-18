function [T,flag] = arReplaceParameterId(T)

global ar
flag = 0;

[~, ia, ib] = intersect(ar.pLabel,cellstr(T.parameterId));
ia2 = setdiff(1:length(ar.pLabel),ia);
ib2 = setdiff(1:size(T.parameterId,1),ib);
pl = ar.pLabel(ia2);
pid = T.parameterId(ib2,:);

for i=1:size(pid,1)
    if strncmp(pid(i,:),'sd',2)
        idx = contains(pl,'noise');
        if any(idx)
            C = strsplit(ar.pLabel{ia2(idx(1))},'_');
            rep = C{1};
        end
        T.parameterId(ib2(i),:) = regexprep(pid(i,:),'sd',rep);
        flag = 1;
    end
end


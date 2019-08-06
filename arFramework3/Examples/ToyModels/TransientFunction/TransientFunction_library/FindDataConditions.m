% Sucht transformationen z.B. epo_level  0.1

function [pold,fp] = FindDataConditions
global ar

pold = cell(size(ar.model));
fp = cell(size(ar.model));
for m=1:length(ar.model)
    for d=1:length(ar.model(m).data)
        pold{m}{d} = ar.model(m).data(d).pold;
        fp{m}{d} = ar.model(m).data(d).fp;
    end
end


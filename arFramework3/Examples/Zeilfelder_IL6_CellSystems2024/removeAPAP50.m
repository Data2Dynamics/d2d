% This dose showed some toxicity and is therefore removed
ds = arFindData('_', 'input_apap', 50);

for a = 1 : numel( ds )
    ar.model(1).data(ds(a)).yExp = NaN;
    ar.model(1).data(ds(a)).yExpStd = NaN;
end

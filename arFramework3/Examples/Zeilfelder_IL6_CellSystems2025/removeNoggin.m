d1 = arFindData('_', 'input_dcf', 500, 'input_noggin', 500);
d2 = arFindData('_', 'input_apap', 10, 'input_noggin', 500);
ds = union(d1, d2);

for a = 1 : numel(ds)
    ar.model.data(ds(a)).yExp = ar.model.data(ds(a)).yExp * NaN;
end

fprintf('Don''t forget to run arLink!\n');
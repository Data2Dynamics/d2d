function arRecompile(varargin)

global ar

ms = {ar.model.name};
ds = cell(size(ms));
for m=1:length(ar.model)
    ds{m} = {ar.model(m).data.name};
    [uni,ia,ib]= unique(regexprep(ds{m},'_nExpID(\d)+',''));
    ds{m} = uni(ib);  % replace zurueck
    ds{m} = ds{m}(sort(ia)); % nur die unique, aber in alter reihenfolge
end

arInit
for m=1:length(ms)
    arLoadModel(ms{m});
end

for m=1:length(ds)
    for d=1:length(ds{m})
%         arLoadData_withoutNormalization(ds{m}{d}, 1,[],[],[]);
        arLoadData(ds{m}{d},m,'xls', 1);
    end
end

arCompileAll(varargin{:})
arSave('Recompile');


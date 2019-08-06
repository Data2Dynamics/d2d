function arInitNEB
% Initializes NEB path method (searches parameters which are used for NEB method)

global ar

ar.merger.neb.steps = length(ar.model)-2; % all except start and end

% find NEB parameters
ar.merger.neb.qFitlable = {};
index = find(contains(ar.pLabel,'_NEB_001'));
for i = 1:length(index)
    mystring = ar.pLabel{index(i)};
    newStr = erase(mystring,'_NEB_001');
    ar.merger.neb.qFitlable{i} = newStr;
end

[~,qFitBasisModel] = ismember(ar.merger.neb.qFitlable,ar.psLabel);
ar.mergerlhs.ps = ar.ps_sorted(:,qFitBasisModel);


end
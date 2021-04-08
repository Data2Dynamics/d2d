% ar = arImportConfig(arIn)
%
% Import ar.config from arIn to global ar

function ar = arImportConfig(arIn)
global ar

fields = fieldnames(arIn.config);
for i = 1:length(fields)
    ar.config.(fields{i}) = arIn.config.(fields{i});
end
arFprintf(1, 'All config options assigned.\n');
end
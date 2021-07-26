function yamlContent = arReadPEtabYaml(name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~contains(name,'.yaml')
    name = [name,'.yaml'];
end

yamlContent = arReadYaml(name);

fns = fieldnames(yamlContent);
for ifn = 1:length(fns)
    if isstr(yamlContent.(fns{ifn}))
        yamlContent.(fns{ifn}) = strrep(yamlContent.(fns{ifn}),'.tsv','');
    end
end


end


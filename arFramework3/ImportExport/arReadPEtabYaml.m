function yamlContent = arReadPEtabYaml(name)
% Reads PEtab yaml file, outputs struct of cell arrays with size 1
% Not supported: more than one tsv file per category, more than one sbml
% file, multiple problems

if ~contains(name,'.yaml')
    name = [name,'.yaml'];
end

%yamlContent = arReadYaml(name);
yamlraw = ReadYaml(name);

if isfield(yamlraw,'problems')
    if length(yamlraw.problems) == 1
        yamlContent.sbml_files = yamlraw.problems{1}.sbml_files;
        yamlContent.condition_files = yamlraw.problems{1}.condition_files;
        yamlContent.measurement_files = yamlraw.problems{1}.measurement_files;
        yamlContent.observable_files = yamlraw.problems{1}.observable_files;
        if isfield(yamlraw,'parameter_file')
            yamlContent.parameter_file = {yamlraw.parameter_file};
        else
            yamlContent.parameter_file = yamlraw.problems{1}.parameter_file;
        end
    else
        error('Yaml files with more than one problem not supported')
    end
else
    yamlContent.sbml_files = yamlraw.sbml_files;
    yamlContent.condition_files = yamlraw.condition_files;
    yamlContent.measurement_files = yamlraw.measurement_files;
    yamlContent.observable_files = yamlraw.observable_files;
    yamlContent.parameter_file = {yamlraw.parameter_file};
end
end


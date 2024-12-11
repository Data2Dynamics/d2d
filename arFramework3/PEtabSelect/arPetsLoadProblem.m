function arPetsLoadProblem(PetsProblemFile)

arguments
    PetsProblemFile (1,:) char = ''
end

% initialize PEtabSelect struct
global pets

% find the problem file (if not provided)
if isempty(PetsProblemFile) || strcmp(PetsProblemFile, "")
    fprintf('Search PEtabSelect folder for .yaml file\n')
    yamlDir = dir(fullfile(pwd(), 'PEtabSelect', '*.yaml'));
    if isempty(yamlDir)
        error('Did not find .yaml file')
    end
    if length(yamlDir)>1
        error('Found more than one .yaml file')
    end
    PetsProblemFile = fullfile(yamlDir.folder, yamlDir.name);
end

% validate that *.yaml file is a valid PEtabSelect problem
validatePetsYaml(PetsProblemFile);

% load the problem file
pets.problem = pets.module.Problem.from_yaml(PetsProblemFile);


function validPetsYaml = validatePetsYaml(filePath)

validPetsYaml = true;

% read yaml file
try
    yamlData = ReadYaml(filePath);
catch
    warning('Could not read PEtabSelect yaml file')
    validPetsYaml = false;
end

% check if yamlData is a 1x1 struct
if ~isstruct(yamlData) || ~all(size(yamlData) == [1,1])
    warning('PEtabSelect yaml file does not contain a single struct')
    validPetsYaml = false;
end

% check if yamlData has the required fields
requiredFields = {'format_version', 'criterion', 'method', 'model_space_files'};
for i = 1:length(requiredFields)
    if ~isfield(yamlData, requiredFields{i})
        warning(['PEtabSelect yaml file does not contain field: ', requiredFields{i}])
        validPetsYaml = false;
    end
end
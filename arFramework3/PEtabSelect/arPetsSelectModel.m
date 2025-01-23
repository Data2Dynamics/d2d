function arPetsSelectModel(PetsProblemFile, options)
% ARPETSSELECTMODEL Run model selection with PEtab Select and d2d.
%
% INPUT:
%   PetsProblemFile (char): Path to the PEtab Select problem file.
%
% NAME-VALUE ARGUMENTS (OPTIONAL):
%   PetsOutputDir (char): Path to the Petab Select output directory.
%   CalibrationLogsDir (char): Path to the d2d fitting logs.
%   d2dFitFunction (function_handle): Function to fit the models in d2d.
%
% OUTPUT:
%   None
%
% SIDE EFFECTS:
%   - A new global variable "pets" (shorthand for "PEtab Select") is
%     created. It contains the python module, the problem object and a
%     history struct.
%   - A folder (specified by name-value argument "PetsOutputDir") is created.
%     It contains logs and results of the model selection.
%   - The ar structs of all calibrated models are saved.
%   - Logs of the d2d calibration are saved in the "CalibrationLogs" folder.
%
% REQUIREMENTS:
%   PEtab Select is a python module called with Matlab's python environment.
%   Instruction for configuration can be found in the d2d wiki (see below).
%
%   The model selection problem must be specified with correspondibg PEtab
%   Select files.
%
%   The candidate model subspaces must either be specified as valid PEtab
%   problems, or alternatively as ar structs (saved to "Results" folder via
%   arSave) with 'dummy' PEtab files.
%
% USAGE:
%   This is the main function of the PEtab Select module. It can be called
%   directly, without any additional initialization.
%
% DOCUMENTATION:
%   - Examples and test cases: arFramework3/Examples/ToyModels/PEtabSelect
%   - d2d wiki: https://github.com/Data2Dynamics/d2d/wiki/Model-selection-with-PEtab-Select
%   - PEtab Select documentation: https://petab-select.readthedocs.io/en/stable/


% required argument
arguments
    PetsProblemFile char ; 
end

% name-value arguments
arguments
    options.PetsOutputDir char = fullfile(pwd(), 'PEtabSelect', 'Results');
    options.CalibrationLogsDir char = fullfile(pwd(), 'CalibrationLogs')
    options.d2dFitFunction function_handle = @() arPetsDefaultFit();
end

%% INITIALIZATION AND IMPORTS
% initialize PEtab Select Module and load the problem
global pets
arPetsInitModule();
loadPetsProblem(PetsProblemFile);

% import some petab_select constants (define shorthands)
CANDIDATE_SPACE = pets.module.constants.CANDIDATE_SPACE;
UNCALIBRATED_MODELS = pets.module.constants.UNCALIBRATED_MODELS;
PREDECESSOR_MODEL = pets.module.constants.PREDECESSOR_MODEL;
TERMINATE = pets.module.constants.TERMINATE;

% import some petab_select functions (define shorthands)
start_iteration = pets.module.ui.start_iteration;
end_iteration = pets.module.ui.end_iteration;
get_best = pets.module.ui.get_best;
models_to_yaml_list = pets.module.models_to_yaml_list;

% initialize empty objects for the first iteration
candidate_space = py.None;
iteration_best_model = py.None;
allTime_best_model = py.None;

%% MODEL SELECTION LOOP
iter = 0;
while iter >= 0
    
    % initialize iteration
    iter = iter + 1;
    iteration = start_iteration(pets.problem, candidate_space);

    % calibrate the predecessor model (if necessary)
    iteration = calibratePredecessor( ...
        iteration, ...
        options.d2dFitFunction, ...
        options.CalibrationLogsDir);
    
    % memorize useful models (as starting point for fitting routine)
    pets.history.current = struct( ...
        'predecessor', getModelHash(iteration{PREDECESSOR_MODEL}), ...
        'iteration_best', getModelHash(iteration_best_model), ...
        'allTime_best', getModelHash(allTime_best_model));

    % calibrate the candidate models
    for idx = 1:length(iteration{UNCALIBRATED_MODELS})
        model = iteration{UNCALIBRATED_MODELS}{idx};
        model = calibrate( ...
            model, ...
            options.d2dFitFunction, ...
            options.CalibrationLogsDir);
        iteration{UNCALIBRATED_MODELS}{idx} = model;
    end

    % finalize iteration
    iteration_results = end_iteration( ...
        pets.problem, ...
        iteration{CANDIDATE_SPACE}, ...
        iteration{UNCALIBRATED_MODELS});
    
    % prepare next iteration
    candidate_space = iteration_results{CANDIDATE_SPACE};

    % check if the model selection should be terminated
    if iteration_results{TERMINATE}
        break
    end

    % get best models (from this iteration and allTime)
    iteration_best_model = get_best(pets.problem, candidate_space.models);
    allTime_best_model = get_best(pets.problem, pets.problem.state.models);
end

%% FINALIZE AND SAVE RESULTS
% create output directory
[~] = mkdir(options.PetsOutputDir);

% save the PEtab Select struct to file
arPetsSave('petsSave', options.PetsOutputDir);

% get the overall best model and save it to yaml file
final_model = get_best(pets.problem, pets.problem.state.models);
final_model.to_yaml(fullfile(options.PetsOutputDir, 'selected_model.yaml'));

% write all calibrated models to yaml file
models_to_yaml_list(pets.problem.state.models, ...
    fullfile(options.PetsOutputDir, 'all_calibrated_models.yaml'));

fprintf('Model selection finished successfully. Results are stored in: \n%s\n', ...
    options.PetsOutputDir)

end


%% SUBFUNCTIONS FOR LOADING PETAB SELECT PROBLEM

% -------------------------------------------------------------------------
function loadPetsProblem(PetsProblemFile)
% Load the PEtab Select problem from a PEtab Select yaml file.

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

end


% -------------------------------------------------------------------------
function validatePetsYaml(filePath)
% Check if PEtab Select yaml file is valid.

% read yaml file
try
    yamlData = ReadYaml(filePath);
catch
    error('Could not read PEtabSelect yaml file.')
end

% check if yamlData is a 1x1 struct
if ~isstruct(yamlData) || ~all(size(yamlData) == [1,1])
    error('PEtabSelect yaml file does not contain a single struct.')
end

% check if yamlData has the required fields
requiredFields = {'format_version', 'criterion', 'method', 'model_space_files'};
missingFields = setdiff(requiredFields, fieldnames(yamlData));
if ~isempty(missingFields)
    error(['PEtabSelect yaml file misses required fields: ', ...
        strjoin(missingFields, ', '), '.'])
end
end


%% SUBFUNCTIONS FOR MODEL CALIBRATION

% -------------------------------------------------------------------------
function iteration = calibratePredecessor(iteration, d2dFitFunction, logfolder)
% Calibrate the predecessor model of the candidates in the iteration.
% This is useful for the arPetsFitWithSmartInitials function.

arguments
    iteration
    d2dFitFunction function_handle
    logfolder char
end

global pets

% get the predecessor model and its hash
PREDECESSOR_MODEL = pets.module.constants.PREDECESSOR_MODEL;
predecessor = iteration{PREDECESSOR_MODEL};
predecessorHash = getModelHash(predecessor);

% check if predecessor model is valid (not empty or virtual)
VIRTUAL_MODEL_HASH = char(py.str(pets.module.model.VIRTUAL_INITIAL_MODEL.hash));
EMTPY_MODEL_HASH = '';
isValidPredecessor = ~any(strcmp(predecessorHash, ...
    {EMTPY_MODEL_HASH, VIRTUAL_MODEL_HASH}));

% check if predecessor model is still uncalibrated
existCalibratedModels = ~isempty(fieldnames(pets.history.calibrated));
isUncalibrated = (~existCalibratedModels || ...
    ~any(strcmp(predecessorHash, {pets.history.calibrated.modelHash})));

% calibrate the predecessor model if possible and necessary
if isValidPredecessor && isUncalibrated
    predecessor = calibrate(predecessor, d2dFitFunction, logfolder);
    iteration{PREDECESSOR_MODEL} = predecessor;
end
end


% -------------------------------------------------------------------------
function model = calibrate(model, d2dFitFunction, logfolder)
% This function calibrates a single PEtab model object.
%
% If a previous fit is available, the model is not fitted again.
% Otherwise, the model translated to an ar struct and calibrated
% using the specified d2dFitFunction.

arguments
    model
    d2dFitFunction function_handle
    logfolder char
end

global ar
global pets

% define the name of the fitted model
modelHash = getModelHash(model);
saveName = sprintf('%s__fitted', modelHash);

% calibration log
[~, ~, ~] = mkdir(logfolder);
diary(fullfile(logfolder, sprintf('%s.log', modelHash)));

try
    % try to load previously fitted model
    arLoad(saveName);
    fprintf('Model %s has been calibrated before. Skip calibration.\n', ...
        modelHash)
catch
    % no previous fit is available -> import and fit the model
    fprintf([ ...
        'Model %s has not been calibrated before. ' ...
        'Import and calibrate it now.\n'], modelHash)
    model2arStruct(model);
    if sum(ar.qFit)~=0
        d2dFitFunction();
    end
    % save the cailbrated d2d workspace
    arSave(saveName, [], 0);
end

% return the fit results to the model object and memorize the model
model = setSelectionCriteria(model);
model = setEstimatedParameters(model);
addModel2History(model);
diary('off');
end


% -------------------------------------------------------------------------
function model2arStruct(model)
% Translate the PEtab Select model object to an ar struct.
%
% If the PEtab problem of the model is already imported and compiled,
% the corresponding ar struct is loaded from the 'Results' folder.
% Otherwise, the PEtab problem is imported with arImportPEtab and
% the ar struct is saved to the 'Results' folder for later use.
%
% The parameters in ar are set according to the PEtab Select model.

global ar
global pets

% try to load d2d version of yaml file
modelSubspace = char(model.model_subspace_id);
yamlFile = char(model.model_subspace_petab_yaml.resolve());
yamlName = getModelSubspacePetabName(yamlFile);
try
    % try to load a workspace from the 'Results' folder
    arLoad(yamlName);
    fprintf('Model subspace %s is available in d2d. Skip Petab Import.\n', ...
        modelSubspace)
catch
    % if not imported and compiled yet, do so now and save for later
    fprintf([ ... 
        'Model subspace %s is not available in d2d. ' ...
        'Import it now from PEtab file: %s.\n'], modelSubspace, yamlFile)
    arImportPEtab(yamlFile);
    arSave(yamlName, [], 0);
end

% set modelHash variable (can be used in fit function or post-processing)
ar.info.petsModelHash = getModelHash(model);

% translate the model parameters to ar struct
ESTIMATE = pets.module.constants.ESTIMATE;
paramStruct = struct(model.parameters);
pLabels = fieldnames(paramStruct);
for ip = 1:length(pLabels)
    pLabel = pLabels{ip};
    value = paramStruct.(pLabel);
    if value == ESTIMATE
        arSetPars(pLabel, [], 1);
    else
        value = double(value);
        idx = arFindPar(pLabel, 'exact');
        if ar.qLog10(idx)
            % the parameters in PEtabSelect are always on linear scale
            % tranform them to log10 scale if necessary
            value = log10(value);
        end
        arSetPars(pLabel, value, 0);
    end
end
end


% -------------------------------------------------------------------------
function model = setEstimatedParameters(model)
% Set the estimated parameters in the PEtab Select model object

global ar

% get logical indices for fitting and log10 transformation
qFit = logical(ar.qFit);
qFitLog10 = logical(ar.qLog10(qFit));

% get labels and values (on a linear scale)
pLabels = ar.pLabel(qFit);
pValues = ar.p(qFit);
pValues(qFitLog10) = 10.^pValues(qFitLog10);

% create a dictionary from the labels and values
paramStruct = cell2struct(num2cell(pValues), pLabels, 2);
estimated_parameters = py.dict(paramStruct);
% set the estimated parameters in the model
model.set_estimated_parameters(estimated_parameters);

end


% -------------------------------------------------------------------------
function model = setSelectionCriteria(model)
% Calculate and set the criteria in the PEtab Select model object
%
% The criteria are calculated in d2d by arGetMerit.
% Alternatively, they could be calculated in petab_select with the
% compute_criterion method of the model.
%
% The forumlas for the criteria AICc and BIC differ slightly betwen PEtab
% Select and d2d. Specifically, the number of priors appears in the PEtab
% Select versions but not in d2d. If we replace "ar.ndata" in the d2d
% expression by "n_measurements + n_priors" from the PEtab Select
% expression, the formulas seem to be identical.
%
% I don't know, if there is an actual difference or which version is
% correct. For now, we use the d2d version.

global ar
global pets

% Calculate the information criteria
arCalcMerit();
[~, merits] = arGetMerit();

% Set criteria in model
criterion = pets.module.constants.Criterion;
model.set_criterion(criterion.AIC,  merits.aic);
model.set_criterion(criterion.AICC, merits.aicc);
model.set_criterion(criterion.BIC,  merits.bic);
model.set_criterion(criterion.NLLH, merits.loglik/2);
model.set_criterion(criterion.SSR,  merits.chi2_res);

end


% -------------------------------------------------------------------------
function addModel2History(model)
% Store some information about the calibrated model in pets.history struct
% This is useful for the arPetsFitWithSmartInitials function.

global ar
global pets

% create model hsitory struct
modelHistory = struct();
modelHistory.modelHash = getModelHash(model);

% get the parameters from the ar struct
modelHistory.pLabel = ar.pLabel;
modelHistory.p = ar.p;
modelHistory.qFit = ar.qFit;
modelHistory.qLog10 = ar.qLog10;

% set the criterion value
CRITERION = pets.problem.criterion;
modelHistory.(char(CRITERION)) = model.get_criterion(CRITERION);

% append to pets history
if isempty(fieldnames(pets.history.calibrated))
    pets.history.calibrated = modelHistory;
else
    pets.history.calibrated(end+1) = modelHistory;
end

end


% -------------------------------------------------------------------------
function yamlName = getModelSubspacePetabName(yamlFile)
% Assign a unique name to the PEtab file of the model subspace.
% This name is stored in pets.history, and is used to save and load the
% workspace to and from the 'Results' folder.

arguments
    yamlFile char
end

global pets

% check if yamlFile has been read before
qYaml = strcmp(yamlFile, pets.history.PetabYamlFiles);
if any(qYaml)
    yamlName = pets.history.PetabNames{qYaml};
    return
end

% define a uniqe name for the yamlFile
% add a number suffix if necessary
[~, yamlName] = fileparts(yamlFile);
uniqueName = yamlName;
suffix = 1;
while ismember(uniqueName, pets.history.PetabNames)
    uniqueName = sprintf('%s_%d', yamlName, suffix);
    suffix = suffix + 1;
end
yamlName = uniqueName;

% add the yamlFile and yamlName to the pets history
pets.history.PetabYamlFiles{end+1} = yamlFile;
pets.history.PetabNames{end+1} = yamlName;
end


% -------------------------------------------------------------------------
function hash = getModelHash(model)
% get the model hash as MATLAB char
if model == py.None
    hash = '';
else
    hash = char(py.str(model.hash));
end
end


% -------------------------------------------------------------------------
function arPetsDefaultFit()
% ARPETSDEFAULTFIT The default fitting routine for model selection.
%
% The routine combines local search with default initial parameter values,
% local serach with initials based on previously calibrated models, and
% a small multistart search with 10 random initial parameters values.

global ar

% Deactivate Bessel Correction
ar.config.useFitErrorCorrection = 0;

% perform fits
arFit();
arPetsFitWithSmartInitials();
arFitLHS(10);
end
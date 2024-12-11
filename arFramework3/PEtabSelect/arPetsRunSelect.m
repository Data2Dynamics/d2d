function arPetsRunSelect(PetsProblemFile, d2dFitFunction)
% run model selection until petab_select terminates

arguments
    PetsProblemFile char = ''
    d2dFitFunction function_handle = @() arPetsDefaultFit()
end

% initialize PEtab Select Module and load the problem
global pets
arPetsInitModule();
arPetsLoadProblem(PetsProblemFile);

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


% initialize an empty candidate space and best-performing models
candidate_space = py.None;
iteration_best_model = py.None;
allTime_best_model = py.None;

% initialize history struct
pets.history = struct();

iter = 0;
while iter >= 0
    
    % initialize iteration
    iter = iter + 1;
    iteration = start_iteration(problem=pets.problem, candidate_space=candidate_space);

    % calibrate the predecessor model (if necessary)
    iteration = maybeCalibratePredecessor(iteration, d2dFitFunction);
    
    % memorize useful models (as starting point for fitting routine)
    pets.history.predecessor = getModelHash(iteration{PREDECESSOR_MODEL});
    pets.history.iteration_best = getModelHash(iteration_best_model);
    pets.history.allTime_best = getModelHash(allTime_best_model);

    % calibrate the candidate models
    for idx = 1:length(iteration{UNCALIBRATED_MODELS})
        model = iteration{UNCALIBRATED_MODELS}{idx};
        model = calibrate(model, d2dFitFunction);
        iteration{UNCALIBRATED_MODELS}{idx} = model;
    end

    % finalize iteration
    iteration_results = end_iteration( ...
        candidate_space=iteration{CANDIDATE_SPACE}, ...
        calibrated_models=iteration{UNCALIBRATED_MODELS});
    
    % prepare next iteration
    candidate_space = iteration_results{CANDIDATE_SPACE};

    % check if the model selection should be terminated
    if iteration_results{TERMINATE}
        break
    end

    % get best models (from this iteration and allTime)
    iteration_best_model = get_best(problem=pets.problem, ...
        models=candidate_space.models);
    allTime_best_model = get_best(problem=pets.problem, ...
        models=candidate_space.calibrated_models.values());

    % todo: print_model(local_best_model)

end

% get the overall best model and save it to yaml file
final_model = get_best(problem=pets.problem, ...
        models=candidate_space.calibrated_models.values());
final_model.to_yaml("selected_model.yaml");

% write all calibrated models to yaml file
models_to_yaml_list(candidate_space.calibrated_models, 'calibratedModels.yaml')

end


function model = calibrate(model, d2dFitFunction, tryLoadingPrevious)
% This function returns a calibrated model object.
%
% If a previous fit is available, it is loaded and returned.
% Otherwise, the model is calibrated using the specified d2dFitFunction.

arguments
    model
    d2dFitFunction function_handle
    tryLoadingPrevious logical = true
end

global ar
global pets

% define the name of the fitted model
modelHash = getModelHash(model);
saveName = sprintf('%s__fitted', modelHash);

% calibration log
logfolder = 'CalibrationLogs';
[~, ~, ~] = mkdir(logfolder);
diary(fullfile(logfolder, sprintf('%s.log', modelHash)));

% try to load previously fitted model
loadSuccess = false;
if tryLoadingPrevious
    loadSuccess = arLoadLatest(saveName);
end

if ~loadSuccess
    % if no previous fit is available, import and fit the model
    model2arStruct(model);
    if sum(ar.qFit)~=0
        d2dFitFunction();
    end
    % save the ar struct
    arSave(saveName);
end

% memorize the current model has and its parameters
addModel2History(modelHash);

% return the fit results to the model object
model = setCriteria(model);
model = setEstimatedParameters(model);

diary('off');

end


function model = setEstimatedParameters(model)
% Set the estimated parameters in the petab_select model object

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


function model = setCriteria(model)
% Calculate the information criteria from ar struct and set them in the model
%
% The criteria are calculated in d2d by arGetMerit.
% Alternatively, they could be calculated in petab_select with the 
% compute_criterion method of the model.
%
% When comparing the codes of d2d and petab_select, I noticed that the forumlas
% for the criteria differ slightly. Specifically, the number of priors appears
% in petab_select versions of BIC and AICc, but not in d2d.
% The formulas are identical if we replace "ar.ndata" by "n_measurements + n_priors".
%
% I don't know, if there is an actual difference. And if so, which version is correct.
% For now, I will use the d2d version, because it is already implemented in arGetMerit.

global ar
global pets

% Calculate the information criteria
arCalcMerit();
[~, merits] = arGetMerit();

% Set criteria in model
criterion = pets.module.constants.Criterion;
model.set_criterion(criterion.AIC,  merits.aic);
model.set_criterion(criterion.AICC, merits.aic);
model.set_criterion(criterion.BIC,  merits.bic);
model.set_criterion(criterion.NLLH, merits.loglik/2);
model.set_criterion(criterion.SSR,  merits.chi2_res);

end


function addModel2History(modelHash)
% Log the parameters of the current model

global ar
global pets

% log the parameters of the current model
modelHistory = struct();
modelHistory.modelHash = modelHash;
modelHistory.pLabel = ar.pLabel;
modelHistory.p = ar.p;
modelHistory.qFit = ar.qFit;
modelHistory.qLog10 = ar.qLog10;

% append to pets history
if ~isfield(pets.history, 'calibrated')
    pets.history.calibrated = modelHistory;
else
    pets.history.calibrated(end+1) = modelHistory;
end

end


function modelStruct = model2arStruct(model)

% initialize d2d
global ar

% translate the model properties to matlab file formats
modelStruct = model2Matlab(model);

% try to load d2d version of yaml file
yamlFile = modelStruct.petab_yaml;
[~, yamlName] = fileparts(yamlFile);
loadSuccess = arLoadLatest(yamlName);

% if not importet and compiled yet, do so now and save for later
if ~loadSuccess
    arImportPEtab(yamlFile);
    saveName = sprintf('%s__import', yamlName);
    arSave(saveName);
end

% set modelHash variable (can be used in fit function or post-processing)
ar.info.petsModelHash = getModelHash(model);

% set parameters correctly
for ip = 1:length(modelStruct.pLabel)

    pLabel = modelStruct.pLabel(ip);
    qFit = modelStruct.qFit(ip);
    qp = strcmp(ar.pLabel, pLabel);

    if sum(qp) == 1
        % setting the value / or qFit based on PEtabSelect Files

        if qFit
            % only activate the fitting, do not change the parameter value
            arSetPars(pLabel, [], 1);
        else
            % set the value and deactivate fitting
            if ar.qLog10(qp)
                % the parameters in PEtabSelect are always on linear scale
                % tranform them to log10 scale if necessary
                arSetPars(pLabel, log10(modelStruct.p(ip)), 0)
            else
                arSetPars(pLabel, modelStruct.p(ip), 0)
            end
        end
        
    elseif sum(qp) == 0
        warning("PEtab_Select parameter '%s' not found in ar struct!", pLabel)
    else
        warning("Multiple d2d parameters with the same label '%s'!", pLabel)
    end
end

end


function modelStruct = model2Matlab(model)
modelStruct = struct(model.to_dict());

% Convert string fields to MATLAB char arrays
% modelStruct.model_id = char(modelStruct.model_id);
% modelStruct.model_subspace_id = char(modelStruct.model_subspace_id);
% modelStruct.model_hash = char(modelStruct.model_hash);
% modelStruct.predecessor_model_hash = char(modelStruct.predecessor_model_hash);
modelStruct.petab_yaml = char(modelStruct.petab_yaml);

% Convert lists to MATLAB arrays
% modelStruct.model_subspace_indices = double(modelStruct.model_subspace_indices);

% Convert parameters to MATLAB and construct ar.p and ar.qFit (subsets)
modelStruct.parameters = struct(modelStruct.parameters);
modelStruct.pLabel = fieldnames(modelStruct.parameters)';
modelStruct.p = NaN(size(modelStruct.pLabel));
modelStruct.qFit = false(size(modelStruct.pLabel));

fields = fieldnames(modelStruct.parameters);
for ip = 1:numel(fields)
    fieldName = fields{ip};
    value = modelStruct.parameters.(fieldName);

    % Directly convert numeric values to matlab double
    if isa(value, 'py.int') || isa(value, 'py.float') || isnumeric(value)
        modelStruct.parameters.(fieldName) = double(value);
        modelStruct.p(ip) = double(value);

    % Convert strings
    elseif isa(value, 'py.str')
        charValue = char(value);
        numValue = str2double(charValue);
        if strcmp(charValue, 'estimate')
            modelStruct.parameters.(fieldName) = charValue;
            modelStruct.qFit(ip) = true;
        elseif ~isnan(numValue)
            modelStruct.parameters.(fieldName) = numValue;
            modelStruct.p(ip) = numValue;
        else
            warning("Unexpected value for parameter '%s'.", fieldName);
        end
    
    else
        typeChar = char(py.getattr(py.type(value), '__name__'));
        warning("Unexpected type %s for parameter '%s'.", typeChar, fieldName);
    end
end

% convert estimated_parameters dictionary
% modelStruct.estimated_parameters = struct(modelStruct.estimated_parameters);
% fields = fieldnames(modelStruct.estimated_parameters);
% for ip = 1:numel(fields)
%     fieldName = fields{ip};
%     value = modelStruct.estimated_parameters.(fieldName);

%     % Directly convert numeric values to matlab double
%     if isa(value, 'py.int') || isa(value, 'py.float')
%         modelStruct.estimated_parameters.(fieldName) = double(value);

%     % Convert strings
%     elseif isa(value, 'py.str')
%         charValue = char(value);
%         numValue = str2double(charValue);
%         if ~isnan(numValue)
%             modelStruct.estimated_parameters.(fieldName) = numValue;
%         else
%             warning("Unexpected value for estimated_parameter '%s'.", fieldName);
%         end
    
%     else
%         typeChar = char(py.getattr(py.type(value), '__name__'));
%         warning("Unexpected type %s for estimated_parameter '%s'.", typeChar, fieldName);
%     end
% end

% modelStruct.criteria = py.dict(modelStruct.criteria);
% modelStruct.criteria = struct(modelStruct.criteria);
    
end



function hash = getModelHash(model)

global pets
VIRTUAL_INITIAL_MODEL = pets.module.constants.VIRTUAL_INITIAL_MODEL;

if model == py.None
    % no model specified
    hash = '';
elseif model == VIRTUAL_INITIAL_MODEL
    % virtual initial model
    hash = char(model);
else
    % get the hash of the model (usual case)
    hash = char(py.str(model.get_hash));
end

end


function iteration = maybeCalibratePredecessor(iteration, d2dFitFunction)
% Calibrate the predecessor model if valid and not already calibrated

global pets

PREDECESSOR_MODEL = pets.module.constants.PREDECESSOR_MODEL;
VIRTUAL_INITIAL_MODEL = pets.module.constants.VIRTUAL_INITIAL_MODEL;

predecessor = iteration{PREDECESSOR_MODEL};
predecessorHash = getModelHash(predecessor);

isValidPredecessor = ~any(strcmp(predecessorHash, {'', char(VIRTUAL_INITIAL_MODEL)}));
isValidHistory = ~isempty(fieldnames(pets.history));
isAlreadyCalibrated = (isValidHistory && any(strcmp(predecessorHash, {pets.history.calibrated.modelHash})));

if isValidPredecessor && ~isAlreadyCalibrated
    predecessor = calibrate(predecessor, d2dFitFunction);
    iteration{PREDECESSOR_MODEL} = predecessor;
end

end


function arPetsDefaultFit()

global ar

ar.config.useFitErrorCorrection = 0;

arFit();
arPetsFitSmartly();
arFitLHS(10);

end
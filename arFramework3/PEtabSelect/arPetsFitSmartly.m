function arPetsFitSmartly(fitPredecessor, fitAllTimeBest, fitIterationBest)
% ARPETSFITSMARTLY Local fits of current ar struct with inital parameters from previous models.
%
% During model selection with arPetsRunSelect, information about calibrated models
% is stored in the global struct "pets.history".
% This information is used to provide initial parameters for local fits of the current model.
% Only the parameters that are estimated in the current model are replaced.
%
% The input arguments can be used to specify which models should be used for the smart fitting.
% By default, all models are used.
%
% Usage: Include this function in the d2dFitFunction provided to arPetsRunSelect.
%        An example is provided in the arPetsDefaultFitFunction.

arguments
    fitPredecessor (1,1) logical = true
    fitAllTimeBest (1,1) logical = true
    fitIterationBest (1,1) logical = true
end

global ar
global pets

if ~isfield(pets, 'history') || isempty(fieldnames(pets.history))
    return
end

% define petab model for virtual initial model
VIRTUAL_INITIAL_MODEL = pets.module.constants.VIRTUAL_INITIAL_MODEL;

% identify useful models (as starting point for fitting routine)
previousModels = {};
if fitPredecessor
    previousModels = [previousModels, pets.history.predecessor];
end
if fitAllTimeBest
    previousModels = [previousModels, pets.history.allTime_best];
end
if fitIterationBest
    previousModels = [previousModels, pets.history.iteration_best];
end

% remove invalid and duplicate entries
previousModels = setdiff(previousModels, {'', char(VIRTUAL_INITIAL_MODEL)});
previousModels = unique(previousModels);

nModels = length(previousModels);
if nModels == 0
    fprintf('No previous models found for smart fitting.\n');
    return
end

% check if calibrated models are available
if ~isfield(pets.history, 'calibrated') || isempty(pets.history.calibrated)
    fprintf('No calibrated models found for smart fitting.\n');
    return
end

% which parameters are fitted
estimatedParameters = ar.pLabel(logical(ar.qFit));

% get the parameters of these models
historyHashes = {pets.history.calibrated.modelHash};
for iModel = 1:nModels
    
    % check if the calibrated model is in history
    idxHistory = strcmp(historyHashes, previousModels{iModel});
    if sum(idxHistory) == 0
        fprintf('Model %s not found in history.\n', previousModels{iModel});
        continue
    elseif sum(idxHistory) > 1
        fprintf('Model %s found multiple times in history.\n', previousModels{iModel});
        continue
    end

    % memorize parameters of current model
    paramsOriginal = ar.p;

    % get parameter struct of previous model
    paramStruct = pets.history.calibrated(idxHistory);

    % replace all estimated parameters with values from previous model
    for iEstim = 1:length(estimatedParameters)
        pLabel = estimatedParameters{iEstim};

        % find parameter in parameter struct of previous model
        qParam = strcmp(pLabel, paramStruct.pLabel);
        if sum(qParam) == 0
            fprintf('Parameter %s not found in model %s.\n', pLabel, previousModels{iModel});
            continue
        elseif sum(qParam) > 1
            fprintf('Parameter %s found multiple times in model %s.\n', pLabel, previousModels{iModel});
            continue
        end

        % set parameter value
        qLog10 = ar.qLog10(arFindPar(pLabel, 'exact'));
        if qLog10 == paramStruct.qLog10(qParam)
            arSetPars(pLabel, paramStruct.p(qParam));
        elseif qLog10 == 1
            arSetPars(pLabel, log10(paramStruct.p(qParam)));
        else % qLog10 == 0
            arSetPars(pLabel, 10^(paramStruct.p(qParam)));
        end

    end

    % check if parameters changed
    if any(ar.p ~= paramsOriginal)
        % perform a local fit
        arFit();
    else
        % previous model did not provide useful initials
        % -> skip fitting
    end

end

end
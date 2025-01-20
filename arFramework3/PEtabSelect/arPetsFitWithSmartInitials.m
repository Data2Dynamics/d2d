function arPetsFitWithSmartInitials(fitPredecessor, fitAllTimeBest, fitIterationBest)
% ARPETSFITWITHSMARTINITIALS Perform local fits with inital parameters from previous models.
%
% INPUT:
%   fitPredecessor   (logical, optional):   Use the predecessor model.
%   fitAllTimeBest   (logical, optional):   Use the all-time best model.
%   fitIterationBest (logical, optional):   Use the iteration best model.
%
% OUTPUT:
%   None
%
% SIDE EFFECTS:
%   The model parameters (ar.p) are modified and optimized with arFit.
%
% EXPLANATION:
%   During model selection with arPetsSelectModel, the estimated parameters
%   of pervious models are stored in the global struct "pets.history".
%   These parameter values are assumed to be good inital values for local
%   fits of the current model. Therefore, the current parameters to be 
%   estimated are initialized with values of previously calibrated models.
%   Then, a local fit is performed with arFit.
%
%   The previous models that are assumed to be beneficial are:
%   - The predecessor model of the current model.
%   - The all-time best model.
%   - The best model of the previous iteration.
%   The user can specify which of these models should be used.
%   Usually, some of these models coincide and only are only used once.
%
% USAGE:
%   It is recommended to include this function in the d2dFitFunction
%   that is provided to arPetsSelectModel. An example is the default  fitting
%   routing for PEtab Select: arPetsSelectModel>arPetsDefaultFitFunction.

arguments
    fitPredecessor (1,1) logical = true
    fitAllTimeBest (1,1) logical = true
    fitIterationBest (1,1) logical = true
end

global ar
global pets

% verify that the history is already initilized and has entries
if isempty(fieldnames(pets.history.calibrated))
    fprintf('No previously calibrated models in pets.histoty, yet.\n');
    return
end

% collect the hashes of the models requested for parameter initilization
previousModels = {};
if fitPredecessor
    previousModels = [previousModels, pets.history.current.predecessor];
end
if fitAllTimeBest
    previousModels = [previousModels, pets.history.current.allTime_best];
end
if fitIterationBest
    previousModels = [previousModels, pets.history.current.iteration_best];
end

% remove invalid and duplicate model hashes
VIRTUAL_MODEL_HASH = char(py.str(pets.module.model.VIRTUAL_INITIAL_MODEL.hash));
EMTPY_MODEL_HASH = '';
previousModels = setdiff(previousModels, {VIRTUAL_MODEL_HASH, EMTPY_MODEL_HASH});
previousModels = unique(previousModels);

% return if valid models exist
nModels = length(previousModels);
if nModels == 0
    fprintf('No previous models found for smart fitting.\n');
    return
end

% use the models from history as initials for parameters
% that should be estimated
estimatedParameters = ar.pLabel(logical(ar.qFit));
historyHashes = {pets.history.calibrated.modelHash};
for iModel = 1:nModels
    
    % check if the calibrated model is in history
    idxHistory = strcmp(historyHashes, previousModels{iModel});
    if sum(idxHistory) == 0
        fprintf('Model %s not found in history.\n', ...
            previousModels{iModel});
        continue
    elseif sum(idxHistory) > 1
        fprintf('Model %s found multiple times in history.\n', ... 
            previousModels{iModel});
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
            fprintf('Parameter %s not found in model %s.\n', ...
                pLabel, previousModels{iModel});
            continue
        elseif sum(qParam) > 1
            fprintf('Parameter %s found multiple times in model %s.\n', ...
                pLabel, previousModels{iModel});
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
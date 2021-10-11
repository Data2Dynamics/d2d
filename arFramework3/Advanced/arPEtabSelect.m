%% Work in Progress
% Optional extension in different function: Write selection problem yaml with input
% Currently only one criterion supported

function arPEtabSelect(venvActPath, yaml, method, limit, initialModel, CalibYamlOut, estimationRoutine, iterationCounter)
if ~exist('iterationCounter') || isempty(iterationCounter)
    iterationCounter = 1;
end
%% Parse function input

% standard settings
if ~exist('selectionProblem') || isempty(yaml)
    yaml = 'selection_problem.yaml';
end
if ~exist('method') || isempty(method)
    method = 'brute_force';
end
if ~exist('level') || isempty(limit)
    limit = 3;
end
if ~exist('initialModel') || isempty(initialModel)
    initialModel = '';
end
if ~exist('CalibYamlOut') || isempty(CalibYamlOut)
    CalibYamlOut = 'calibrated_it_001.yaml';
end
if ~exist('estimationRoutine') || isempty(estimationRoutine)
    estimationRoutine = @arFit;
end

if exist('venvActPath') && ~isempty(venvActPath)
    initstr = sprintf('source %s; ', venvActPath);
else
    initstr = '';
    venvActPath = '';
end

%% Check if petab_select installation
syscom = [initstr, 'petab_select --help'];
[status,~] = system(syscom);
if status ~= 0
    error(sprintf('Calling petab_select from the command line failed.\nPlease check your Python environment and the PEtab-select installation.'))
end

%% Call PEtab-select to generate candidate models
fprintf('arPEtabSelect: Generating candidate models...\n')
SelectionProblem = ReadYaml(yaml);

syscom = [initstr,...
    'petab_select candidates ',  ...
    ' -y ', yaml, ...
    ' -s output', filesep, 'state.dill',...
    ' -o output', filesep, 'models.yaml', ...
    ' -m ', method, ...
    ' -l ', num2str(limit), ...
    ];
if ~isempty(initialModel)
    syscom = [syscom, ' -b ' ,initialModel];
end
[status,cmdout] = system(syscom);

if status ~= 0
    error(sprintf('Error while running petab_select from command line.\n Command line message:\n %s',cmdout)); %#ok<SPERR>
end

%% Process candidate models
iterationCounter %debug
CandidateModels = ReadYaml(['output' filesep 'models.yaml']);
nModels = size(CandidateModels,2);

if nModels == 0
    fprintf('arPEtabSelect: Finished after iteration %i - no candidate models found.\n', iterationCounter-1)
    return
end
fprintf('arPEtabSelect: Calibrating candidate models...\n')

for jModel = 1:nModels
    % Load & compile
    arInit
    doPreEq = false;
    arImportPEtab(CandidateModels{jModel}.petab_yaml,doPreEq)
    ar.config.useFitErrorCorrection = 0;
    
    % Import parameter settings
    pars = fieldnames(CandidateModels{jModel}.parameters);
    estimatedPars = {};
    for iPar = 1:length(pars)
        if CandidateModels{jModel}.parameters.(pars{iPar}) == 'estimate'
            arSetPars(pars{iPar},[],1)
            estimatedPars{end+1} = pars{iPar};
        else
            parId = arFindPar(pars{iPar});
            parValue = CandidateModels{jModel}.parameters.(pars{iPar});
            if ar.qLog10(parId) == 1
                arSetPars(pars{iPar},log10(parValue),0)
            else
                arSetPars(pars{iPar},parValue,0)
            end
        end
    end
    
    % Estimate
    %estimationRoutine;

    if sum(ar.qFit) > 0
       arFit
       arFitLHS(10)
    end
    arCalcMerit
    [~, allmerits] = arGetMerit;
    
    % Collect criteria
    switch SelectionProblem.criterion
        case 'AIC'
            theCriterion = allmerits.aic;
        case 'AICc'
            theCriterion = allmerits.aicc;
        case 'BIC'
            theCriterion = allmerits.bic;
        case 'LogL'
            theCriterion = allmerits.loglik;
        otherwise
            error('Invalid criteria for model %i: %s',jModel,CandidateModels{jModel}.petab_yaml)
    end
    
    %Write results
    calibCands{jModel}.criteria.(SelectionProblem.criterion) = theCriterion;
    calibCands{jModel}.model_id = CandidateModels{jModel}.model_id;
    calibCands{jModel}.parameters = CandidateModels{jModel}.parameters;
    calibCands{jModel}.petab_yaml = CandidateModels{jModel}.petab_yaml;

    if isempty(estimatedPars)
        calibCands{jModel}.estimated_parameters = 'null';
    else
        for iPar = 1:length(estimatedPars)
            calibCands{jModel}.estimated_parameters.(estimatedPars{iPar}) = 10^ar.p(arFindPar(estimatedPars{iPar}))*ar.qLog10(arFindPar(estimatedPars{iPar})) + ar.p(arFindPar(estimatedPars{iPar}))*(1-ar.qLog10(arFindPar(estimatedPars{iPar})));
        end
    end
end
WriteYaml(CalibYamlOut,calibCands);

%% Next iteration
fprintf('arPEtabSelect: Iteration %i complete. Continuing with next iteration\n',iterationCounter)
nextIterationYamlOut = sprintf('calibrated_it_%03i.yaml',iterationCounter+1);
arPEtabSelect(venvActPath, yaml, method, limit, CalibYamlOut, nextIterationYamlOut, estimationRoutine,iterationCounter+1)
end
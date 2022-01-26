% arPEtabSelect(venvActPath, yaml, limit, estimRoutine)
%
% Run model selection using PEtab-select (https://github.com/PEtab-dev/petab_select)
% Requires installation of the petab_select Python3 package. arPEtabSelect
% uses the command line interface (CLI).
%
%   [venvActPath]   Path to a python virtual environment (venv) activation
%                   (https://docs.python.org/3/library/venv.html). If left
%                   empty, the CLI will use the global python installation.
%
%   [yaml]          Path to the PEtab-select problem yaml-file 
%                   [petab_select_problem.yaml]
%
%   [limit]         Maximum number of processed candidate models per
%                   iteration. Cannot be used with method 'brute_force'
%                   [3]
%
%   [estimRoutine]  D2D commands for parameter estimation 
%                   ['arFit; arFitLHS(10)']
%
% Examples:
%   Process selection problem in '~/test_cases/0001/petab_select_problem'
%   with a python venv in '~/d2d_python_venv/bin/activate'. Use 0.3 as
%   initial guess for all estimated parameters and additionally perform
%   multi-start optimization with 20 runs:
%
%   arPEtabSelect('~/d2d_python_venv/bin/activate', ...
%       '~/test_cases/0001/petab_select_problem',...
%       3, 
%       'ar.p(ar.qFit == 1) = 0.3; arFit; arFitLHS(20)')
% 
%
% (Leave empty initialModel & iterationCtr. Those are internal arguments
% to allow for recursive function call)

function arPEtabSelect(venvActPath, yaml, limit, estimRoutine, initialModel, iterationCtr)
if ~exist('iterationCtr') || isempty(iterationCtr)
    iterationCtr = 1;
end
if ~exist('petab-select', 'dir')
    mkdir('petab-select')
end

%% Parse function input

% standard settings
if ~exist('selectionProblem') || isempty(yaml)
    yaml = 'petab_select_problem.yaml';
end
if ~exist('initialModel') || isempty(initialModel)
    initialModel = '';
end
if ~exist('estimRoutine') || isempty(estimRoutine)
    estimRoutine = 'arFit; arFitLHS(10);';
end
if exist('venvActPath') && ~isempty(venvActPath)
    initstr = sprintf('source %s; ', venvActPath);
else
    initstr = '';
    venvActPath = '';
end

% if method = brute_force, do not allow limit argument
SelectionProblem = ReadYaml(yaml);
if strcmp(SelectionProblem.method,'brute_force')
    if ~exist('limit') || isempty(limit)
        limit = '';
    else
        error('Limit argument not allowed if method = brute_force')
    end
else
    if ~exist('limit') || isempty(limit)
        limit = 3;
    end
end

terminateFlag = 0;
CalibYamlOut = ['petab-select', filesep, sprintf('calibrated_it_%03i.yaml',iterationCtr)];

%% Check if petab_select installation
syscom = [initstr, 'petab_select --help'];
[status,~] = system(syscom);
if status ~= 0
    error(sprintf('Calling petab_select from the command line failed. ...\nPlease check your Python environment and the PEtab-select installation.'))
end

%% Call PEtab-select to generate candidate models
fprintf('arPEtabSelect: Generating candidate models...\n')
syscom = [initstr,...
    'petab_select candidates ',  ...
    ' -y ', yaml, ...
    ' -s output', filesep, 'state.dill',...
    ' -o output', filesep, 'models.yaml', ...
    ' --relative-paths ', ...
    ];
if ~isempty(initialModel)
    syscom = [syscom, ' -b ', initialModel];
end
if ~isempty(limit)
    syscom = [syscom, ' -l ', num2str(limit)];
end

[status,cmdout] = system(syscom);

if status ~= 0
    error(sprintf('Error while running petab_select candidates from command line.\n Command line message:\n %s',cmdout)); %#ok<SPERR>
end

%% Process candidate models
%iterationCtr %debug
CandidateModels = ReadYaml(['output' filesep 'models.yaml']);
nModels = size(CandidateModels,2);

if nModels < 1
    fprintf('arPEtabSelect: Finished after iteration %i - no (more) candidate models found.\n', iterationCtr-1)
    terminateFlag = 1;
end
if terminateFlag == 0
    fprintf('arPEtabSelect: Calibrating candidate models...\n')
    
    for jModel = 1:nModels
        % Load & compile
        arInit
        doPreEq = false;
        arImportPEtab(['output', filesep, CandidateModels{jModel}.petab_yaml],doPreEq)
        ar.config.useFitErrorCorrection = 0;
        
        % Import parameter settings
        pars = fieldnames(CandidateModels{jModel}.parameters);
        estimatedPars = {};
        
        for iPar = 1:length(pars)
            parIndex(iPar) = find(ismember(ar.pLabel,pars(iPar)));
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
        
        % Add all estimated parameters that are not in model yaml
        if sum(ar.qFit) > 0
            for iModPar = 1:length(ar.qFit)
                % If not already treated above
                if sum(iModPar == parIndex) == 0
                    estimatedPars{end+1} = ar.pLabel{iModPar};
                end
            end
        end
        
        % Estimate        
        if sum(ar.qFit) > 0
            eval(estimRoutine);
        end
        arCalcMerit
        [~, allmerits] = arGetMerit;
        
        % Collect criteria
        criteria.AIC = allmerits.aic;
        criteria.AICc = allmerits.aicc;
        criteria.BIC = allmerits.bic;
        %criteria.nllh = allmerits.loglik/(2);
        
        calibCands{jModel}.criteria = criteria;
        calibCands{jModel}.model_id = CandidateModels{jModel}.model_id;
        calibCands{jModel}.parameters = CandidateModels{jModel}.parameters;
        calibCands{jModel}.petab_yaml = CandidateModels{jModel}.petab_yaml;
        
        calibCands{jModel}.model_subspace_id = CandidateModels{jModel}.model_subspace_id;
        calibCands{jModel}.model_hash = CandidateModels{jModel}.model_hash;
        calibCands{jModel}.predecessor_model_hash = CandidateModels{jModel}.predecessor_model_hash;
        calibCands{jModel}.model_subspace_indices = CandidateModels{jModel}.model_subspace_indices;
        
        if isempty(estimatedPars)
            calibCands{jModel}.estimated_parameters = struct();
        else
            for iPar = 1:length(estimatedPars)
                calibCands{jModel}.estimated_parameters.(estimatedPars{iPar}) = ...
                    10^ar.p(arFindPar(estimatedPars{iPar}))*ar.qLog10(arFindPar(estimatedPars{iPar})) + ...
                    ar.p(arFindPar(estimatedPars{iPar}))*(1-ar.qLog10(arFindPar(estimatedPars{iPar})));
            end
        end
    end
    WriteYaml(CalibYamlOut,calibCands);
    
    %% Find best model of current iteration
    syscom = [initstr,...
        'petab_select best ', ...
        ' -y ', yaml,...
        ' -m ', CalibYamlOut,...
        ' -o petab-select', filesep, sprintf('best_model_it_%03i.yaml',iterationCtr),...
        ' -s output', filesep, 'state.dill',...
        ' --relative-paths ',...
        ];
    [status,cmdout] = system(syscom);
    if status ~= 0
        error(sprintf('Error while running petab_select best from command line.\n Command line message:\n %s',cmdout)); %#ok<SPERR>
    end
    
    %% Read current and previous iteration's criterion (and stop)
    if iterationCtr > 1
        prevIt = ReadYaml(['petab-select',filesep,...
            sprintf('best_model_it_%03i.yaml',iterationCtr-1)]);
        prevItCrit = prevIt.criteria.(SelectionProblem.criterion);
        currentIt = ReadYaml(['petab-select',filesep,...
            sprintf('best_model_it_%03i.yaml',iterationCtr)]);
        currentItCrit = currentIt.criteria.(SelectionProblem.criterion);
        
        if round(currentItCrit,4) > round(prevItCrit,4)
            fprintf('arPEtabSelect: Finished after iteration %i - criterion worse than in iteration %i.\n', iterationCtr, iterationCtr-1)
            terminateFlag = 1;
        end
    end
    
    if strcmp(SelectionProblem.method,'brute_force')
        terminateFlag = 2;
    end
end

%% Write results yaml
if terminateFlag == 1
    copyfile(['petab-select',filesep,...
        sprintf('best_model_it_%03i.yaml',iterationCtr-1)],...
        ['petab-select',filesep,'selected_model.yaml'])
    return
elseif terminateFlag == 2
    copyfile(['petab-select',filesep,...
        sprintf('best_model_it_%03i.yaml',iterationCtr)],...
        ['petab-select',filesep,'selected_model.yaml'])
    return
end

%% Next iteration
fprintf('arPEtabSelect: Iteration %i complete. Continuing with next iteration\n',iterationCtr)
arPEtabSelect(venvActPath, yaml, limit, estimRoutine, CalibYamlOut ,iterationCtr+1)
end
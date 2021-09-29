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

initstr = '';
if exist('venvActPath') && ~isempty(venvActPath)
    initstr = sprintf('source %s; ', venvActPath);
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
    arInit
    doPreEq = false;
    arImportPEtab(CandidateModels{jModel}.petab_yaml,doPreEq)
    estimationRoutine;
    [~, allmerits] = arGetMerit;
            
    % ADAPT FOR MULTIPLE CRITERIA
    switch SelectionProblem.criterion
        case 'AIC'
             theCriterion = allmerits.loglik;  
        case 'LogL'
             theCriterion = allmerits.loglik;  
        otherwise
             error('Invalid criteria for model %i: %s',jModel,CandidateModels{jModel}.petab_yaml)
    end
end

calibCands = CandidateModels;
for jModel = 1:nModels
    calibCands{jModel}.criteria.(SelectionProblem.criterion) = theCriterion;
    calibCands{jModel}.estimated_parameters = 'null';
    % add values of estimated parameters
end 
WriteYaml(CalibYamlOut,calibCands);

%% Next iteration
fprintf('arPEtabSelect: Iteration %i complete. Continuing with next iteration\n',iterationCounter)
nextIterationYamlOut = sprintf('calibrated_it_%03i.yaml',iterationCounter+1);
arPEtabSelect(venvActPath, yaml, method, limit, CalibYamlOut, nextIterationYamlOut, estimationRoutine,iterationCounter+1)
end
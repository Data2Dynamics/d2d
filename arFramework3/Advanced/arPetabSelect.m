%% Work in Progress
% Optional extension in different function: Write selection problem yaml with input
% Currently only one criterion supported

function arPetabSelect(venvActPath)
%% Parse function input





%% Check if petab_select is installed in system commands and give error if not yet installed
initstr = '';
if exist('venvActPath') && ~isempty(venvActPath)
    initstr = sprintf('source %s; ', venvActPath);
end

command = [initstr, 'petab_select --help'];
[status,~] = system(command);
if status ~= 0
    error(sprintf('petab_select could not be called as a system command.\nPlease make sure to install petab_select using pip install petab_select on your machine.'))
end



%% Run petab_select code with system commands
SelectionProblem = ReadYaml('selection_problem.yaml');

% Create output folder? Currently still missing
if ~isfolder('output')
    mkdir('output')
end
command = [initstr,...
'petab_select candidates ',  ...
' -y selection_problem.yaml ', ...
' -s output', filesep, 'state.dill',...
' -o output', filesep, 'models.yaml', ...
' -m brute_force', ...
' -l 3'
];
[status,cmdout] = system(command);

if status ~= 0
    error(sprintf('Error while running petab_select from command line.\n Command line message:\n %s',cmdout)); %#ok<SPERR>
end

%% Import models.yaml written by petab_select script
CandidateModels = ReadYaml(['output' filesep 'models.yaml']);


%% Run d2d fits & calculate criterion values with d2d functions
nModels = size(CandidateModels,2);
for jModel = 1:nModels
    % collect PEtab files 
    bb = ReadYaml(CandidateModels{jModel}.petab_yaml);
    %lastfilesepPos = find(CandidateModels(jModel).petab_yaml == filesep, 1, 'last');
    %PeTsvPath = CandidateModels(jModel).petab_yaml(1:lastfilesepPos);
    %PeTsvPath = cellfun(@(x) [PeTsvPath,x],{bb.sbml_files, bb.observable_files, bb.measurement_files, bb.condition_files, bb.parameter_file},'UniformOutput',false);
    
    arInit
    doPreEq = false;
    
    arImportPEtab(bb,doPreEq)
    % Fit model 
    arFit
    [~, allmerits] = arGetMerit;
    
    % Save criteria
    % ADAPT FOR MULTIPLE CRITERIA
    switch SelectionProblem.criterion
        case 'AICc'
             bb(jModel).criteria = allmerits.loglik;  
        case 'LogL'
             bb(jModel).criteria = allmerits.loglik;  
        otherwise
             error('Invalid criteria for model nr%i:\n%s',jModel,CandidateModels(jModel).petab_yaml)
    end
    
    
    
end

%% Write calibrated first iteration .yaml including criterion values




%% use petab_select candidates -b calibrated first iteration.yaml with method (forward selection...)





%% Final output 

end
%% Setup the d2d models and save them

% define model.def files and names for the two d2d models
modelsDefs = {'model1', 'model2'};
modelNames = {'d2dModel1', 'd2dModel2'};
petabFiles = [];

for isub = 1:2
    % create the d2d model
    arInit;
    arLoadModel(modelsDefs{isub});
    arLoadData('data');
    arCompileAll;

    % save the workspace and PEtab Files with matching names
    % -> This is important if you want to avoid overwriting the original d2d
    % problem with a version that is first exported to PEtab and then reimported.
    arSave(modelNames{isub}, [], 0);
    arExportPEtab(modelNames{isub})
    petabFiles(end+1) = string(fullfile(pwd(), 'PEtab', sprintf('%s.yaml', modelNames{isub})));
end

%% Create the PETab Select problem files 

% setup directory and file names
petsFilesDir = fullfile(pwd(), 'PEtabSelect');
petsOutputDir = fullfile(petsFilesDir, 'output');
[~] = mkdir(petsFilesDir);
[~] = mkdir(petsOutputDir);
petsProblemFile = fullfile(petsFilesDir, 'petab_select_problem.yaml');
petsModelSpaceFile = fullfile(petsFilesDir, 'model_space.tsv');

% create problem .yaml file
petsProblem = struct();
petsProblem.format_version = 'beta_1';
petsProblem.criterion = 'AIC';
petsProblem.method = 'forward';
petsProblem.model_space_files = {petsModelSpaceFile};
WriteYaml(petsProblemFile, petsProblem);

% create model space .tsv table
petsModelSpace = table();

% model subspace M1 (no annihilations of x1 and x2, encoded in model1.def)
petsModelSpace.model_subspace_id(1) = "M1";
petsModelSpace.model_subspace_petab_yaml(1) = petabFiles(1);
petsModelSpace.(ar.pLabel{1})(1) = "0;0.2;estimate"; 
petsModelSpace.(ar.pLabel{2})(1) = "0;0.1;estimate";

% model subspace M2 (estimate annihilation rate k3, encoded in model2.def)
petsModelSpace.model_subspace_id(2) = "M2";
petsModelSpace.model_subspace_petab_yaml(2) = petabFiles(2);
petsModelSpace.(ar.pLabel{1})(2) = "0;0.2;estimate";
petsModelSpace.(ar.pLabel{2})(2) = "0;0.1;estimate";

writetable(petsModelSpace, petsModelSpaceFile,...
    'Delimiter', '\t', 'FileType', 'text')


%% Run PEtabSelect
arPetsSelectModel(petsProblemFile, PetsOutputDir=petsOutputDir);
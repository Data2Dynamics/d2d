%% Setup the d2d model and save it
arInit;
arLoadModel('model');
arLoadData('data');
arCompileAll;

% save the workspace and PEtab Files with matching names
% -> This is important if you want to avoid overwriting the original d2d
% problem with a version that is first exported to PEtab and then reimported.
modelName = 'd2dModel';
arSave(modelName, [], 0);
arExportPEtab(modelName);
petabFile = string(fullfile(pwd(), 'PEtab', sprintf('%s.yaml', modelName)));


%% Create the PEtab Select problem files

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
% add single model subspace "M"
petsModelSpace.model_subspace_id(1) = "M";
petsModelSpace.model_subspace_petab_yaml(1) = petabFile;
petsModelSpace.(ar.pLabel{1})(1) = "0;0.2;estimate"; 
petsModelSpace.(ar.pLabel{2})(1) = "0;0.1;estimate"; 
petsModelSpace.(ar.pLabel{3})(1) = "0;estimate";
writetable(petsModelSpace, petsModelSpaceFile,...
    'Delimiter', '\t', 'FileType', 'text')


%% Run PEtabSelect
arPetsSelectModel(petsProblemFile, PetsOutputDir=petsOutputDir);
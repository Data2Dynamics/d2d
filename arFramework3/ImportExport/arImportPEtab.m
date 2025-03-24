
%%
% ToDo
% - Ergänzungen: events, precalib

function arImportPEtab(PEtabName, doPreEq, model, dataFolder, dataFilenames)
% arImportPEtabNew(PEtabName, transformData, doPreEq, dataFilenames, model)
%   Imports PEtab files and data files into d2d arFramework
%   PEtabName:      Name of PEtab file (without extension)
%   doPreEq:        Perform pre-equilibration [default: true]
%   model:          Name of d2d model file [default: create from PEtab xml file]
%   dataFolder:     Folder where data files are stored [default: DataPEtab]
%   dataFilenames:  Names of data files to be imported [default: create from PEtab measurements file]

global ar

if(isempty(ar))
    fprintf('No ar struct found. Initializing d2d.\n')
    arInit();
end
if isfield(ar, 'model')
    fprintf('ar struct already contains a model. Re-initializing d2d.\n')
    arInit();
end

if ~exist('PEtabName','var') || isempty(PEtabName)
    % find .yaml files in folder PEtab
    folder = 'PEtab';
    files = dir(fullfile(folder, '*.yaml'));
    % count number of yaml files
    if isempty(files)
        error('Did not find .yaml file')
    elseif length(files)>1
        error('Found more than one .yaml file')
    end
    PEtabName = files(1).name;
end

if ~exist('doPreEq','var') || isempty(doPreEq)
    doPreEq = true;
end

if ischar(PEtabName)
    % import from yaml file
    yamlDir = dir(['PEtab' filesep strrep(PEtabName,'.yaml','') '.yaml']);
    if isempty(yamlDir) % exists .yaml file?
        yamlDir = dir(fullfile('PEtab', [strrep(PEtabName,'.yaml','') '.yaml']));
        if isempty(yamlDir) % exists .yaml file?
            error('Did not find .yaml file')
        end
    end
    PEtabName = fullfile('PEtab', yamlDir.name);
    yamlContent = arReadPEtabYaml(PEtabName);
    yamlPath = yamlDir.folder;
    fprintf('Import .yaml file with name: %s\n',PEtabName);

    % check number of files per category, only one per category allowed atm
    petabFiles = {'sbml_files', 'observable_files', 'measurement_files',...
        'condition_files', 'parameter_file'};
    inputArgs = {};
    for i = 1:numel(petabFiles)
        [out,numberOfEls] = extractFromStruct(yamlContent, petabFiles{i});
        if numberOfEls > 1
            error('arImportPEtab: YAML file can only refer to one file per category.')
        end
        inputArgs{end+1} = out;
    end
    if ~exist('model','var')
        model = [];
    end

    %%
    if (~exist('dataFolder','var') || isempty(dataFolder)) & (~exist('dataFilenames','var') || isempty(dataFilenames))
        dataFolder = 'DataPEtab';
        if ~exist(dataFolder, 'dir')  % Check if folder does not exist
            mkdir(dataFolder);  % Create folder
        end
        dataFilenames = arWriteDataFromPEtab(PEtabName, dataFolder);
    elseif (exist('dataFolder','var') || ~isempty(dataFolder)) & (~exist('dataFilenames','var') || isempty(dataFilenames))
        files = dir(fullfile(dataFolder, '*.xls'));
        dataFilenames = arrayfun(@(f) f.name(1:end-4), files, 'UniformOutput', false);
    elseif (exist('dataFolder','var') || ~isempty(dataFolder)) & (exist('dataFilenames','var') || ~isempty(dataFilenames))
        files = dir(fullfile(dataFolder, dataFilenames,[strrep(PEtabName,'.xls','') '.xls']));
        dataFilenames = arrayfun(@(f) f.name(1:end-4), files, 'UniformOutput', false);
    end

    %%
    arImportPEtab(cellfun(@(x) [yamlPath, filesep, x], [inputArgs{:}], ...
        'UniformOutput', false), doPreEq, model, dataFolder, dataFilenames);
    return
end

%%
PEobs = [strrep(PEtabName{2},'.tsv',''),'.tsv'];
PEmeas = [strrep(PEtabName{3},'.tsv',''),'.tsv'];
PEconds = [strrep(PEtabName{4},'.tsv',''),'.tsv'];
PEparas = [strrep(PEtabName{5},'.tsv',''),'.tsv'];

if ~exist('modelname','var') || isempty(modelname)
    sbmlmodel = dir(PEtabName{1});
    if isempty(sbmlmodel)
        error('No SBML file found!');
    else
        ar.petab.sbml = sbmlmodel;
    end
    if length(sbmlmodel)>1
        error('Found more than one SBML model file');
    end
    [~,modelname,~] = arParseSBML([sbmlmodel.folder filesep sbmlmodel.name]);
end

%% ar
arLoadModel(modelname);
splitedFilenames = split(dataFilenames,filesep);
for i=1:length(dataFilenames)
    arLoadData(char(splitedFilenames(1,i,2)),[],[],[],'DataPath',char(splitedFilenames(1,i,1)));
end
arCompileAll;
ar.config.fiterrors = 1;
arLoadParsPEtab(PEparas);
arFindInputs(); % might overwrite parameters due to ar.pExtern, but input times might be in parameters table.
arLoadParsPEtab(PEparas);

%%
%% pre-equilibration
% if doPreEq
%     tstart = -1e7;
%     for imodel = 1:length(ar.model)
%         Tdat = T{imodel,1};
%
%         if isfield(table2struct(Tdat), 'preequilibrationConditionId')
%             uniqueSimConds = unique(Tdat.simulationConditionId);
%             % uniqueSimConds =
%             uniquePreEqConds = unique(Tdat.preequilibrationConditionId);
%             % uniquePreEqConds =
%             if isa(uniquePreEqConds, 'string')
%                if all(strcmp(uniquePreEqConds, ""))
%                    uniquePreEqConds = [];
%                end
%             else
%                 if all(isnan(uniquePreEqConds))
%                     uniquePreEqConds = [];
%                 end
%             end
%
%             if numel(uniquePreEqConds) > 1
%                 error('More than one pre-equiblibration condition currently not supported.')
%             end
%
%             for iPreEqCond = 1:size(uniquePreEqConds,1)
%                 preEqCondId = convertStringsToChars(uniquePreEqConds(iPreEqCond));
%                 preEqCond = arFindCondition(preEqCondId, 'conservative');
%                 simConds = [];
%                 for iSimCond = 1:size(uniqueSimConds,1)
%                     simCondId = convertStringsToChars(uniqueSimConds(iSimCond));
%                     simCond = arFindCondition(simCondId, 'conservative');
%                     simCondDat = Tdat(Tdat.simulationConditionId==simCondId, :);
%                     if all(simCondDat.preequilibrationConditionId==preEqCondId)
%                         simConds = [simConds, simCond];
%                     end
%                     % simConds(end+1) = arFindCondition(convertStringsToChars(uniqueSimConds(iSimCond)), 'conservative');
%                     % addSimConds = find(cellfun(@(x) ~strcmp(x, convertStringsToChars(uniquePreEqConds(iPreEqCond))), {ar.model.data.name}));
%                     % simConds = [simConds, addSimConds];
%                 end
%                 simConds = unique(simConds);
%                 % simConds
%                 arSteadyState(imodel, preEqCond, simConds, tstart);
%             end
%             for iSimCond = 1:length(uniqueSimConds)
%                 Tcondi = Tcond(Tcond.conditionId == uniqueSimConds(iSimCond),:);
%
%                 % iSimCond auf richtige provozieren
%                 %iSimCondAr = ar.model(m).data(iSimCond).cLink;
%                 iSimCondAr = find(cellfun(@(x) strcmp(x, uniqueSimConds(iSimCond)), {ar.model(m).data.name}));
%                 iSimCondAr = ar.model(m).data(iSimCondAr).cLink;
%                 % alt:
%                 % iSimCondAr = find(cellfun(@(x) strcmp(x, uniqueSimConds(iSimCond)), {ar.model(m).data.name}));
%
%                 for iCol = 1:length(Tcondi.Properties.VariableNames)
%                     idxState = find(strcmp(Tcondi.Properties.VariableNames{iCol},ar.model.xNames));
%                     if ~isempty(idxState)
%                         arAddEvent(m,iSimCondAr,0.0001,ar.model.x{idxState}, 0, Tcondi.(ar.model.xNames{idxState}));
%                         % ar = arAddEvent([ar], model, condition, timepoints, [statename], [A], [B],  [sA], [sB])
%                     end
%                 end
%             end
%         end
%     end
% end
%
% % events
% for iev = 1:length(eventStruct)
%     arAddEvent(1, 'all', eventStruct(iev).time, ...clce
%         eventStruct(iev).state, 0, eventStruct(iev).value);
% end
% end

end

function [out,numberOfEls] = extractFromStruct(struct, field)
if iscell(struct.(field))
    out = struct.(field);
    numberOfEls = numel(out);
else
    out = {struct.(field)};
    numberOfEls = 1;
end
end
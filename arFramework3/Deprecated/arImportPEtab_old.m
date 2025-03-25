function arImportPEtab(name, doPreEq)
% arImportPEtab(name, doPreEq)
% Import parameter estimation problem formulated in the PEtab standard.
%
%   name     Path to PEtab yaml file
%
%            Cell array of paths to model, observables, measurements,
%            conditions and parameters file (in this order):
%               arImportPEtab({'mymodel', 'myobs', 'mymeas', 'mycond', 'mypars'})
%
%            If empty, use *.yaml file from PEtab folder (if exactly one exists)
%
%   doPreEq  Apply pre-equilibration if specified in PEtab files [true]
%
% See also
%       arExportPEtab
%
% References
%   - PEtab standard: https://petab.readthedocs.io/en/latest/
%   - D2D Wiki: https://github.com/Data2Dynamics/d2d/wiki/Support-for-PEtab
global ar

if(isempty(ar))
    fprintf('No ar struct found. Initializing d2d.\n')
    arInit();
end
if isfield(ar, 'model')
    fprintf('ar struct already contains a model. Re-initializing d2d.\n')
    arInit();
end
if ~exist('doPreEq','var') || isempty(doPreEq)
    doPreEq = true;
end

%TODO: read in multiple sbmls & save these paths in ar struct

%% Import from yaml file
if ~exist('name', 'var') || isempty(name)
    % find yaml file in PEtab folder
    fprintf('Search PEtab folder for .yaml file\n')
    yamlDir = dir(fullfile('PEtab', '*.yaml'));
    if isempty(yamlDir)
        error('Did not find .yaml file')
    end
    if length(yamlDir)>1
        error('Found more than one .yaml file')
    end
    name = fullfile('PEtab', yamlDir.name);
end
if ischar(name)
    % import from yaml file 
    yamlDir = dir([strrep(name,'.yaml','') '.yaml']);
    if isempty(yamlDir) % exists .yaml file?
        yamlDir = dir(fullfile('PEtab', [strrep(name,'.yaml','') '.yaml']));
        if isempty(yamlDir) % exists .yaml file?
            error('Did not find .yaml file')
        end
    end
    name = fullfile('PEtab', yamlDir.name);
    yamlContent = arReadPEtabYaml(name);
    yamlPath = yamlDir.folder;
    fprintf('Import .yaml file with name: %s\n',name);
    
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
    
    arImportPEtab(cellfun(@(x) [yamlPath, filesep, x], [inputArgs{:}], 'UniformOutput', false),doPreEq);
    % also check arReadPEtabYaml
    return
end

%% Import from list of PEtab files
sbmlmodel = dir([strrep(name{1},'.xml','') '.xml']);
if isempty(sbmlmodel)
    error('No SBML file found!');
else
    ar.petab.sbml = sbmlmodel;
end
if length(sbmlmodel)>1
    error('Found more than one SBML model file')
    %out = stringListChooser({sbmlmodel.name});
    %sbmlmodel= sbmlmodel(out);% set sbmlmodel to chosen
end

[~,modelname,eventStruct] = arParseSBML([sbmlmodel.folder filesep sbmlmodel.name]);
arLoadModel(modelname); % arLoadModel(strrep(sbmlmodel.name,'.xml','')); % arParseSBML has flag overwrite with defaule false

ar.config.backup_modelAndData = false;

PEobs = [strrep(name{2},'.tsv',''),'.tsv'];
PEmeas = [strrep(name{3},'.tsv',''),'.tsv'];
PEconds = [strrep(name{4},'.tsv',''),'.tsv'];
PEparas = [strrep(name{5},'.tsv',''),'.tsv'];

T = cell(2, length(ar.model));
for m = 1:length(ar.model)
    [T{m,1}, T{m,2}] = arLoadDataPEtab(PEmeas,PEobs,m);
end
Tcond = arLoadCondPEtab(PEconds, T);

% Compilation
arCompileAll
ar.config.fiterrors = 1;

arLoadParsPEtab(PEparas);

arFindInputs(); % might overwrite parameters due to ar.pExtern, but input times might be in parameters table.
arLoadParsPEtab(PEparas);

%% pre-equilibration
if doPreEq
    tstart = -1e7;
    for imodel = 1:length(ar.model)
        Tdat = T{imodel,1};
        Tobs = T{imodel,2};
        
        if isfield(table2struct(Tdat), 'preequilibrationConditionId')
            uniqueSimConds = unique(Tdat.simulationConditionId);
            % uniqueSimConds = 
            uniquePreEqConds = unique(Tdat.preequilibrationConditionId);
            % uniquePreEqConds = 
            if isa(uniquePreEqConds, 'string')
               if all(strcmp(uniquePreEqConds, ""))
                   uniquePreEqConds = [];
               end
            else
                if all(isnan(uniquePreEqConds))
                    uniquePreEqConds = [];
                end
            end
            
            if numel(uniquePreEqConds) > 1
                error('More than one pre-equiblibration condition currently not supported.')
            end
            
            for iPreEqCond = 1:size(uniquePreEqConds,1)
                preEqCondId = convertStringsToChars(uniquePreEqConds(iPreEqCond));
                preEqCond = arFindCondition(preEqCondId, 'conservative');
                simConds = [];
                for iSimCond = 1:size(uniqueSimConds,1)
                    simCondId = convertStringsToChars(uniqueSimConds(iSimCond));
                    simCond = arFindCondition(simCondId, 'conservative');
                    simCondDat = Tdat(Tdat.simulationConditionId==simCondId, :);
                    if all(simCondDat.preequilibrationConditionId==preEqCondId)
                        simConds = [simConds, simCond];
                    end
                    % simConds(end+1) = arFindCondition(convertStringsToChars(uniqueSimConds(iSimCond)), 'conservative');
                    % addSimConds = find(cellfun(@(x) ~strcmp(x, convertStringsToChars(uniquePreEqConds(iPreEqCond))), {ar.model.data.name}));
                    % simConds = [simConds, addSimConds];
                end
                simConds = unique(simConds);
                % simConds
                arSteadyState(imodel, preEqCond, simConds, tstart);
            end
            for iSimCond = 1:length(uniqueSimConds)
                Tcondi = Tcond(Tcond.conditionId == uniqueSimConds(iSimCond),:);
                
                % iSimCond auf richtige provozieren
                %iSimCondAr = ar.model(m).data(iSimCond).cLink;
                iSimCondAr = find(cellfun(@(x) strcmp(x, uniqueSimConds(iSimCond)), {ar.model(m).data.name}));
                iSimCondAr = ar.model(m).data(iSimCondAr).cLink;
                % alt:
                % iSimCondAr = find(cellfun(@(x) strcmp(x, uniqueSimConds(iSimCond)), {ar.model(m).data.name}));

                for iCol = 1:length(Tcondi.Properties.VariableNames)
                    idxState = find(strcmp(Tcondi.Properties.VariableNames{iCol},ar.model.xNames));
                    if ~isempty(idxState)
                        arAddEvent(m,iSimCondAr,0.0001,ar.model.x{idxState}, 0, Tcondi.(ar.model.xNames{idxState}));
                        % ar = arAddEvent([ar], model, condition, timepoints, [statename], [A], [B],  [sA], [sB])
                    end
                end
            end
        end
    end
end

% events
for iev = 1:length(eventStruct)
    arAddEvent(1, 'all', eventStruct(iev).time, ...
        eventStruct(iev).state, 0, eventStruct(iev).value);
end
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
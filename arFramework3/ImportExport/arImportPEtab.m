function arImportPEtab(name, doPreEq)
% arImportPEtab(name, doPreEq)
% Import parameter estimation problem formulated in the PEtab standard.
%
%   name     Path to PEtab yaml file
%
%            Alternatively: Cell array of paths to model, observables, measurements,
%            conditions and parameters file (in this order):
%               arImportPEtab({'mymodel', 'myobs', 'mymeas', 'mycond', 'mypars'})
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
    error('Please initialize by arInit')
end
if isfield(ar, 'model')
    error('Please initialize by arInit')
end
if ~exist('doPreEquilibration','var') || isempty(doPreEq)
    doPreEq = true;
end

%TODO: read in multiple sbmls & save these paths in ar struct
if ischar(name) %yaml
    yamlDir = dir([strrep(name,'.yaml','') '.yaml']);
    if isempty(yamlDir) % exists .yaml file?
        error('Did not find yaml file')
    end
    yamlContent = arReadPEtabYaml(name);
    yamlPath = yamlDir.folder;
    fprintf('Found yaml file with name %s\n',name)
    
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
    
    arImportPEtab(cellfun(@(x) [yamlPath, filesep, x], [inputArgs{:}], 'UniformOutput', false),doPreEq)
    % also check arReadPEtabYaml
    return
else % no yaml file
    sbmlmodel = dir([strrep(name{1},'.xml','') '.xml']);
end

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

arParseSBML([sbmlmodel.folder filesep sbmlmodel.name])
arLoadModel(strrep(sbmlmodel.name,'.xml',''))

PEobs = [strrep(name{2},'.tsv',''),'.tsv'];
PEmeas = [strrep(name{3},'.tsv',''),'.tsv'];
PEconds = [strrep(name{4},'.tsv',''),'.tsv'];
PEparas = [strrep(name{5},'.tsv',''),'.tsv'];

T = cell(2, length(ar.model));
for m = 1:length(ar.model)
    [T{m,1}, T{m,2}] = ...
        arLoadDataPEtab(PEmeas,PEobs,m);
end
Tcond = arLoadCondPEtab(PEconds, T);

% Compilation
arCompileAll
ar.config.fiterrors = 1;

arLoadParsPEtab(PEparas);

arFindInputs % might overwrite parameters due to ar.pExtern, but input times might be in parameters table.
arLoadParsPEtab(PEparas);

if exist('Tms','var')
    [BothPars,ia] = intersect(ar.pLabel,Tms_fn);
    for i=1:length(BothPars)
        if strcmp(Tms.(BothPars{i})(ms),'-')       % if Tms.(BothPars{i})(ms) == 0
            arSetPars(BothPars(i), ar.p(ia), 0);
            %else
            %    arSetPars(BothPars(i), ar.p(ia)*Tms.(BothPars{i})(ms));
        end
    end
end

% pre-equilibration
if doPreEq
    tstart = -1e7;
    for imodel = 1:length(ar.model)
        Tdat = T{imodel,1};
        Tobs = T{imodel,2};
        
        if isfield(table2struct(Tdat), 'preequilibrationConditionId')
            uniqueSimConds = unique(Tdat.simulationConditionId);
            uniquePreEqConds = unique(Tdat.preequilibrationConditionId);
            if ~strcmp(class(uniquePreEqConds),'string')
                if all(isnan(uniquePreEqConds))
                    uniquePreEqConds = [];
                end
            end
            
            if numel(uniquePreEqConds) > 1
                error('More than one pre-equiblibration condition currently not supported.')
            end
            
            for ipreeqcond = 1:size(uniquePreEqConds,1)
                preEqCond = arFindCondition(convertStringsToChars(uniquePreEqConds(ipreeqcond)), 'conservative');
                simConds = [];
                for isimcond = 1:size(uniqueSimConds,1)
                    %simConds(end+1) = arFindCondition(convertStringsToChars(uniqueSimConds(isimcond)), 'conservative');
                    simConds(end+1) = ...
                        find(cellfun(@(x) ~strcmp(x, convertStringsToChars(uniquePreEqConds(ipreeqcond))), {ar.model.data.name}));
                end
                arSteadyState(imodel, preEqCond, simConds, tstart)
            end
            for isimu = 1:length(uniqueSimConds)
                Tcondi = Tcond(Tcond.conditionId == uniqueSimConds(isimu),:);
                iSimuAr = find(cellfun(@(x) strcmp(x, uniqueSimConds(isimu)), {ar.model(m).data.name}));
                for iCol = 1:length(Tcondi.Properties.VariableNames)
                    idxState = find(strcmp(Tcondi.Properties.VariableNames{iCol},ar.model.xNames));
                    if ~isempty(idxState)
                        arAddEvent(m,iSimuAr,0.0001,ar.model.x{idxState}, 0, Tcondi.(ar.model.xNames{idxState}))
                        % ar = arAddEvent([ar], model, condition, timepoints, [statename], [A], [B],  [sA], [sB])
                    end
                end
            end
        end
    end
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
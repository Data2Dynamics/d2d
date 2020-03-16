function arImportPEtab(name, doPreEq)
% arImportPEtab(name, doPreEq)
% Import parameter estimation problem formulated in the PEtab standard. 
%
%   name     Cell array of filenames for model, observables, measurements,
%            conditions and parameters file (in this order):
%               arImportPEtab({'mymodel', 'myobs', 'mymeas', 'mycond', 'mypars'})
%
%            Alternatively, a string can be provided that together with the 
%            keywords 'model', 'observables', 'measurements', 'conditions', 
%            'parameters' uniquely defines a set of PEtab files by the
%            pattern 'name*[keyword]'
%               arImportPEtab('my')
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

% init PEtab fields




%TODO: read in multiple sbmls & save these paths in ar struct
if ~exist('name','var') || isempty(name) 
    sbmlmodel = dir(['**' filesep '*.xml']);
else
    if ischar(name)
        sbmlmodel = dir(['**' filesep '*' name '*.xml']);
    else
        sbmlmodel = dir(['**' filesep name{1}]);
        if isempty(sbmlmodel)
            sbmlmodel = dir(['**' filesep name{1} '.xml']);
        end
    end
end
if isempty(sbmlmodel)
    error('No SBML file found! Switch your path to working directory.');
end
if length(sbmlmodel)>1
    out = stringListChooser({sbmlmodel.name});
    sbmlmodel= sbmlmodel(out);                 % set sbmlmodel to chosen
end

if ~exist('name','var') || isempty(name) 
    pe_dir = dir('**/*.tsv');
elseif ~ischar(name)
    pe_dir = dir(['**/*' name{2} '*.tsv']);
else
    pe_dir = dir(['**/*' name '*.tsv']);
end

if isempty(sbmlmodel)
    error('No tsv files found! Switch your path to working directory.');
end
pe_dir = pe_dir(1).folder;

arParseSBML([sbmlmodel.folder filesep sbmlmodel.name])
arLoadModel(strrep(sbmlmodel.name,'.xml',''))

% make dir case sensitive!
if exist('name','var') && ~isempty(name) && ischar(name)
        PEobs = dir([pe_dir filesep name '*observables' '.tsv']);
        PEmeas = dir([pe_dir filesep name '*measurements' '.tsv']);
        PEconds = dir([pe_dir filesep name '*conditions' '.tsv']);
        PEparas = dir([pe_dir filesep name '*parameters' '.tsv']);
else
    if ~exist('name','var') || isempty(name) || length(name)<2 || isempty(name{2})
        PEobs = dir([pe_dir filesep sprintf('*%s*.tsv', 'observables')]);    
    else
        PEobs = dir([name{2} '.tsv']);
    end
    if ~exist('name','var') || isempty(name) || length(name)<3 || isempty(name{3})
        PEmeas = dir([pe_dir filesep sprintf('*%s*.tsv', 'measurements')]);
    else
        PEmeas = dir([name{3} '.tsv']);
    end
    if ~exist('name','var') || isempty(name) || length(name)<4 || isempty(name{4})
        PEconds = dir([pe_dir filesep sprintf('*%s*.tsv', 'conditions')]);
    else
        PEconds = dir([name{4} '.tsv']);
    end
    if ~exist('name','var') || isempty(name) || length(name)<5 || isempty(name{5})
        PEparas = dir([pe_dir filesep sprintf('*%s*.tsv', 'parameters')]);
    else
        PEparas = dir([name{5} '.tsv']);
    end
end

if length(PEparas) > 1 || length(PEmeas) > 1 || length(PEconds) > 1 || length(PEobs) > 1
    error('Found more than one peTAB file. Please provide a unique name substring or specify names of the individual files directly.')
end

% ToDo: Loop over several models
T = cell(2, length(ar.model));
for m = 1:length(ar.model)
    [T{m,1}, T{m,2}] = ...
        arLoadDataPEtab([pe_dir filesep PEmeas.name],[pe_dir filesep PEobs.name],m);
end
Tcond = arLoadCondPEtab([pe_dir filesep PEconds.name], T);

% Compilation
arCompileAll
ar.config.fiterrors = 1;

arLoadParsPEtab([pe_dir filesep PEparas.name]); 
arFindInputs % might overwrite parameters due to ar.pExtern, but input times might be in parameters table.
arLoadParsPEtab([pe_dir filesep PEparas.name]); % get para values from parameter label

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
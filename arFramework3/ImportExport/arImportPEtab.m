% arImportPEtab(filename, doPreEquilibration)
%
% Import PEtab model and data format to Data2Dynamics. 
%
%       name                String that must be contained in the filenames
%                           of all files to load
%       doPreEquilibration  Apply pre-equilibration if specified in PEtab files [true]
%
% See also
%       arExportPEtab
%
% References
%   - https://github.com/ICB-DCM/PEtab/blob/master/doc/documentation_data_format.md

function arImportPEtab(name, doPreEquilibration)
global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if isfield(ar, 'model')
    error('Existing ar workspace detected. Please initialize by arInit for PEtab import')
end

if ~exist('name','var') || isempty(name)
    name = '';
else
    name = [name '*'];
end

if ~exist('doPreEquilibration') || isempty(doPreEquilibration)
    doPreEquilibration = true;
end

%TODO: read in multiple sbmls & save these paths in ar struct
sbmlmodel = dir(['**' filesep '*' name '.xml']);
if isempty(sbmlmodel)
    error('No sbml file found! Switch your path to working directory.');
end
if length(sbmlmodel)>1
    out = stringListChooser({sbmlmodel.name});
    sbmlmodel= sbmlmodel(out);                 % set sbmlmodel to chosen
end

pe_dir = dir('**/*.tsv');
if isempty(sbmlmodel)
    error('No tsv files found! Switch your path to working directory.');
end
pe_dir = pe_dir(1).folder;

arParseSBML([sbmlmodel.folder filesep sbmlmodel.name])
arLoadModel(strrep(sbmlmodel.name,'.xml',''))

% make dir case sensitive!
PEobs = dir([pe_dir filesep sprintf('%s%s*.tsv', name, 'obs')]);
PEmeas = dir([pe_dir filesep sprintf('%s%s*.tsv', name, 'meas')]);
PEconds = dir([pe_dir filesep sprintf('%s%s*.tsv', name, 'cond')]);
PEparas = dir([pe_dir filesep sprintf('%s%s*.tsv', name, 'par')]);

if length(PEparas) > 1 || length(PEmeas) > 1 || length(PEconds) > 1 || length(PEobs) > 1
    error('Found more than one peTAB file.')
end

% ToDo: Loop over several models
T = cell(2, length(ar.model));
for m = 1:length(ar.model)
    [T{m,1}, T{m,2}] = ...
        arLoadDataPEtab([pe_dir filesep PEmeas.name],[pe_dir filesep PEobs.name],m);
end
Tcond = arLoadCondPEtab([pe_dir filesep PEconds.name]);


% Compilation
arCompileAll

arLoadParsPEtab([pe_dir filesep PEparas.name]);
arFindInputs

% pre-equilibration
if doPreEquilibration
    for imodel = 1:length(ar.model)
        Tdat = T{imodel,1};
        Tobs = T{imodel,2};
        
        if isfield(Tdat, 'preEquilibrationId')
            uniqueSimConds = unique(Tdat.simulationConditionId);
            uniquePreEqConds = unique(Tdat.preEquilibrationId);
            
            for ipreeqcond = 1:size(uniquePreEqConds,1)
                preEqCond = arFindCondition(convertStringsToChars(uniquePreEqConds(ipreeqcond)));
                simConds = [];
                for isimcond = 1:size(uniqueSimConds,1)
                    simConds(end+1) = arFindCondition(convertStringsToChars(uniqueSimConds(isimcond)));
                end
                arSteadyState(imodel, preEqCond, simConds)
            end
        end
    end
end
end
% arImportPEtab(name, doPreEquilibration, tstart)
%
% Import PEtab model and data format to Data2Dynamics. 
%
%       name                Character string of model name or
%                           cell array of input files               []
%                           'model1' or {'namemodel1.xml','_OBS_model1','_MEAS_model1','_COND_model1','_PARS_model'}
%                           if name is empty, pattern search of *.xml *obs*.tsv *cond*.tsv *meas*.tsv *pars*.tsv is performed
%       doPreEquilibration  Apply pre-equilibration if specified in PEtab files [true]
%       tstart              Starting time for pre-equilibration [0]
%
% See also
%       arExportPEtab
%
% References
%   - https://github.com/ICB-DCM/PEtab/blob/master/doc/documentation_data_format.md
%
% Examples:
% arImportPEtab
% arImportPEtab({'namemodel1_l2v4.xml','_OBS_model1.tsv'})
% arImportPEtab({'namemodel1_l2v4.xml',[],'_PARS_model.tsv'})

function arImportPEtab(name, doPreEquilibration)
global ar

if(isempty(ar))
    error('please initialize by arInit')
end
if isfield(ar, 'model')
    error('Existing ar workspace detected. Please initialize by arInit for PEtab import')
end
if ~exist('doPreEquilibration','var') || isempty(doPreEquilibration)
    doPreEquilibration = true;
end    

%TODO: read in multiple sbmls & save these paths in ar struct
if ~exist('name','var') || isempty(name) 
    sbmlmodel = dir(['**' filesep '*.xml']);
else
    if ischar(name)
        sbmlmodel = dir(['**' filesep '*' name '.xml']);
    else
        sbmlmodel = dir(['**' filesep name{1}]);
        if isempty(sbmlmodel)
            sbmlmodel = dir(['**' filesep name{1} '.xml']);
        end
    end
end
if isempty(sbmlmodel)
    error('No sbml file found! Switch your path to working directory.');
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
        PEobs = dir([pe_dir filesep 'observables_' name '.tsv']);
        PEmeas = dir([pe_dir filesep 'measurementData_' name '.tsv']);
        PEconds = dir([pe_dir filesep 'experimentalCondition_' name '.tsv']);
        PEparas = dir([pe_dir filesep 'parameters_' name '.tsv']);
else
    if ~exist('name','var') || isempty(name) || length(name)<2 || isempty(name{2})
        PEobs = dir([pe_dir filesep sprintf('*%s*.tsv', 'obs')]);    
    else
        PEobs = dir([name{2} '.tsv']);
    end
    if ~exist('name','var') || isempty(name) || length(name)<3 || isempty(name{3})
        PEmeas = dir([pe_dir filesep sprintf('*%s*.tsv', 'meas')]);
    else
        PEmeas = dir([name{3} '.tsv']);
    end
    if ~exist('name','var') || isempty(name) || length(name)<4 || isempty(name{4})
        PEconds = dir([pe_dir filesep sprintf('*%s*.tsv', 'cond')]);
    else
        PEconds = dir([name{4} '.tsv']);
    end
    if ~exist('name','var') || isempty(name) || length(name)<5 || isempty(name{5})
        PEparas = dir([pe_dir filesep sprintf('*%s*.tsv', 'par')]);
    else
        PEparas = dir([name{5} '.tsv']);
    end
end

if length(PEparas) > 1 || length(PEmeas) > 1 || length(PEconds) > 1 || length(PEobs) > 1
    error('Found more than one peTAB file.')
end

% ToDo: Loop over several models
T = cell(2, length(ar.model));
for m = 1:length(ar.model)
    [T{m,1}, T{m,2}] = ...
        arLoadDataPEtab([pe_dir filesep PEmeas.name],[pe_dir filesep PEobs.name],m);
end
arLoadCondPEtab([pe_dir filesep PEconds.name]);

% Compilation
arCompileAll

arLoadParsPEtab([pe_dir filesep PEparas.name]); 
arFindInputs % might overwrite parameters due to ar.pExtern, but input times might be in parameters table.
arLoadParsPEtab([pe_dir filesep PEparas.name]); % get para values from parameter label

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
                arSteadyState(imodel, preEqCond, simConds, tstart)
            end
        end
    end
end
    
end
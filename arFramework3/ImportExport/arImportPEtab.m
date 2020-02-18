% arImportPEtab(filename)
%
% Import Parameter estimation problem to Data2Dyanmics
%
% This is required to load additional information about parameters from
% PEtab data standard
%
% Example
%       arImportPEtab('good_fit')
%
% See also
%       arExportPEtab
%
% References
%   - https://github.com/ICB-DCM/PEtab/blob/master/doc/documentation_data_format.md

function arImportPEtab(name)
global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if isfield(ar, 'model')
    error('Existing ar workspace detected. Please initialize by arInit for PEtab import')
end

if ~exist('name') || isempty(name)
    name = '*';
end

%TODO: read in multiple sbmls & save these paths in ar struct
sbmlmodel = dir(['**/' name '.xml']);
if isempty(sbmlmodel)
    error('No sbml file found! Switch your path to working directory.');
end
out = stringListChooser({sbmlmodel.name});
sbmlmodel= sbmlmodel(out);                 % set sbmlmodel to chosen

pe_dir = dir('**/*.tsv');
if isempty(sbmlmodel)
    error('No tsv files found! Switch your path to working directory.');
end
pe_dir = pe_dir(1).folder;

arParseSBML([sbmlmodel.folder filesep sbmlmodel.name])
arLoadModel(strrep(sbmlmodel.name,'.xml',''))

% make dir case sensitive!
PEobs = dir([pe_dir filesep sprintf('*%s*.tsv', 'obs')]);
PEmeas = dir([pe_dir filesep sprintf('*%s*.tsv', 'meas')]);
PEconds = dir([pe_dir filesep sprintf('*%s*.tsv', 'cond')]);
PEparas = dir([pe_dir filesep sprintf('*%s*.tsv', 'par')]);

if length(PEparas) > 1 || length(PEmeas) > 1 || length(PEconds) > 1 || length(PEobs) > 1
    error('Found more than one peTAB file.')
end

% ToDo: Loop over several models

for m = 1:length(ar.model)
    arLoadDataPEtab([pe_dir filesep PEmeas.name],[pe_dir filesep PEobs.name],m);
end

arLoadCondPEtab([pe_dir filesep PEconds.name]);

% Compilation
arCompileAll

arLoadParsPEtab([pe_dir filesep PEparas.name]);

end
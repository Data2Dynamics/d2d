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
    name = '';
end

%TODO: make this dynamic OR save these paths in ar struct

pe_dir = 'PEtab/';
sbml_dir = 'SBML_singleConditions/';

counter = 0;
sbmlmodel = dir([pe_dir filesep name  '*.xml']);

if length(sbmlmodel) ~= 1
    error('Not exactly one .xml file found.')
end

arParseSBML([sbmlmodel.folder filesep sbmlmodel.name])
arLoadModel(strrep(sbmlmodel.name,'.xml',''))

% make dir case sensitive!
PEobs = dir([pe_dir name sprintf('*%s*.tsv', '_OBS_')]);
PEmeas = dir([pe_dir name sprintf('*%s*.tsv', '_MEAS_')]);
PEconds = dir([pe_dir name sprintf('*%s*.tsv', '_COND_')]);
PEparas = dir([pe_dir name sprintf('*%s*.tsv', '_PARS_')]);


if length(PEparas) > 1 || length(PEmeas) > 1 || length(PEconds) > 1 || length(PEobs) > 1
    error('Found more than one peTAB file.')
end

% ToDo: Loop over several models

for m = 1:length(ar.model)
    arLoadDataPEtab([pe_dir PEmeas.name],[pe_dir PEobs.name],m);
end

arLoadCondPEtab([pe_dir PEconds.name]);

% Compilation
arCompileAll


arLoadParsPEtab([pe_dir PEparas.name]);

end
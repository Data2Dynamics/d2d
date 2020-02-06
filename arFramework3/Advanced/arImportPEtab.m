% arImportPEtab(filename)
%
% Import Parameter estimation problem to Data2Dyanmics
% !! arInit before arImportPEtab !!
%
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

if ~exist('name') || isempty(name)
    name = '';
end

%TODO: make this dynamic OR save these paths in ar struct

pe_dir = 'PEtab/';
sbml_dir = 'SBML_singleConditions/';

counter = 0;
sbmlmodel = dir([pe_dir '*.xml']);
if isempty(sbmlmodel)
    sbmlmodel = dir([pe_dir sbml_dir name '*.xml']);
end 

if length(sbmlmodel) ~= 1
    error('Not exactly one .xml file found.')
end

arImportSBML([sbmlmodel.folder filesep sbmlmodel.name])

arInit
arLoadModel(strrep(sbmlmodel.name,'.xml',''))

% make dir case sensitive!
PEobs = dir([pe_dir name sprintf('*%s*.tsv', '_OBS_')]);
PEmeas = dir([pe_dir name sprintf('*%s*.tsv', '_MEAS_')]);
PEconds = dir([pe_dir name sprintf('*%s*.tsv', '_COND_')]);
PEparas = dir([pe_dir name sprintf('*%s*.tsv', '_PARS_')]);


if length(PEparas) ~= 1
    error('Not exactly one parameter specification file found.')
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
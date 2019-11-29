% arImportPEtab(filename)
%
% Import Parameter estimation problem to Data2Dyanmics
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

arImportSBML([sbmlmodel.folder '/' sbmlmodel.name])

arInit
arLoadModel(strrep(sbmlmodel.name,'.xml',''))

PEdatas = dir([pe_dir name sprintf('*%s*.tsv', 'meas_')]);
PEconds = dir([pe_dir name sprintf('*%s*.tsv', 'cond_')]);
PEparas = dir([pe_dir name sprintf('*%s*.tsv', 'pars_')]);

if length(PEparas) ~= 1
    error('Not exactly one parameter specification file found.')
end

for j = 1:length(PEdatas)
    arLoadDataPEtab([pe_dir PEdatas(j).name]);
end

for j = 1:length(PEconds)
    arLoadCondPEtab([pe_dir PEconds(j).name]);
end

arCompileAll

arLoadParsPEtab([pe_dir PEparas.name]);

end
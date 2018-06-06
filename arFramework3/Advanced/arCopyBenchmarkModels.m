% arCopyBenchmarkModels(subset,target_path)
% 
%   Copies a specific set of example models (benchmarks) to a local folder.
%   
% 
%   subset      Default: '20'
%               The benchmark models published in Hass et al, 2018)
% 
%               '5'
%               The five benchmark models used in Kreutz et al, 2016)
% 
% 
%   target_path     Default: pwd

function arCopyBenchmarkModels(subset,target_path)
if ~exist('subset','var') || isempty(subset)
    subset = '20';
elseif isnumeric(subset)
    subset = num2str(subset);
end

if ~exist('target_path','var') || isempty(target_path)
    target_path = pwd;
end


switch lower(subset)
    case {'20',''} % default
        models = {
            'Bachmann_MSB2011',...
            'Becker_Science2010',...
            'Beer_MolBiosyst2014',...
            'Boehm_JProteomeRes2014',...
            'Brannmark_JBC2010',...
            'Bruno_Carotines_JExpBio2016',...
            'Chen_MSB2009',...
            'Crauste_ImmuneCells_CellSystems2017',...
            'Fiedler_BMC2016',...
            'Fujita_SciSignal2010',...
            'Isensee_JCB2018',...
            'Lucarelli_TGFb_2017',...
            'Merkle_JAK2STAT5_PCB2016',...
            'Raia_CancerResearch2011',...
            'Reelin_PONE2017',...
            'Schwen_InsulinMouseHepatocytes_PlosOne2014',...
            'Sobotta_Frontiers2017',...
            'Swameye_PNAS2003',...
            'Weber_BMC2015',...
            'Zheng_PNAS2012'};
    case '5'
        models = {
            'Bachmann_MSB2011',...
            'Becker_Science2010',...
            'Beer_MolBiosyst2014',...
            'Raia_CancerResearch2011',...
            'Swameye_PNAS2003'};        
    otherwise
        error('Subset unknown');
end

DoCopy(models, target_path, subset);





function DoCopy(models, target_path, subset)

source_path = [fileparts(which('arInit')),filesep,'Examples'];

for i=1:length(models)
    fprintf('Copying %s ...\n',models{i});
    copyfile([source_path,filesep,models{i}],[target_path,filesep,models{i}]);
    
    switch lower(subset) % model-specific operatios which has to be done:
        case {'20',''}
            if strcmp(models{i},'Merkle_JAK2STAT5_PCB2016')==1  % just keep SetupFinal.m
                delete([target_path,filesep,models{i},filesep,'SetupCFUE.m']);
                delete([target_path,filesep,models{i},filesep,'SetupComprehensive.m']);
                delete([target_path,filesep,models{i},filesep,'SetupH838.m']);
                delete([target_path,filesep,models{i},filesep,'SetupSens.m']);
            end
    end
end



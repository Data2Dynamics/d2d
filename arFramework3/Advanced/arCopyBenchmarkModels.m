% arCopyBenchmarkModels([subset], [target_path])
% 
%   Copies a specific set of example models (benchmarks) to a local folder.
% 
%   subset      ['20'] default value corresponding to the
%               Benchmark models published in Hass et al, 2018
% 
%               '5'
%               The five benchmark models used in Kreutz et al, 2016
% 
%               The (starting characters of) model names can also be
%               specified, e.g.:
%               arCopyBenchmarkModels('Becker')
%               arCopyBenchmarkModels('Swa')
% 
%   target_path     Default: pwd
% 
% Examples
% arCopyBenchmarkModels              % 20 default benchmark models
% arCopyBenchmarkModels('Bachmann')  % only the Bachmann model
% arCopyBenchmarkModels('Ra')        % only the Raia model (abbreviation Ra only fits for Raia)
% 
% 
% Reference:  
% Helge Hass, Carolin Loos, Elba Raimundez Alvarez, Jens
% Timmer, Jan Hasenauer, Clemens Kreutz, Benchmark Problems for Dynamic
% Modeling of Intracellular Processes doi: https://doi.org/10.1101/404590  
% 
% C. Kreutz
% New Concepts for Evaluating the Performance of Computational Methods.
% IFAC-PapersOnLine (2016) 49(26): 63-70.

function arCopyBenchmarkModels(subset,target_path)
if ~exist('subset','var') || isempty(subset)
    subset = '20';
elseif isnumeric(subset)
    subset = num2str(subset);
end

if ~exist('target_path','var') || isempty(target_path)
    target_path = pwd;
end

models20 = {
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


switch lower(subset)
    case {'20',''} % default
        models = models20;
    case '5'
        models = {
            'Bachmann_MSB2011',...
            'Becker_Science2010',...
            'Beer_MolBiosyst2014',...
            'Raia_CancerResearch2011',...
            'Swameye_PNAS2003'};        
    otherwise
        ind = strmatch(lower(subset),lower(models20));
        if ~isempty(ind)
            models = models20(ind);
        else
            error('Subset unknown');
        end
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
            if strcmp(models{i},'Lucarelli_TGFb_2017')==1  % just keep reduced model with genes
                rmdir([target_path,filesep,models{i},filesep,'TGFb_ComplexModel_ModelSelection_L1'],'s');
                rmdir([target_path,filesep,models{i},filesep,'TGFb_ComplexModel_withGenes'],'s');
                movefile([target_path,filesep,models{i},filesep,'TGFb_ComplexModel_WithGenes_Reduced',filesep,'*'],...
                    [target_path,filesep,models{i},filesep,'.']);
                rmdir([target_path,filesep,models{i},filesep,'TGFb_ComplexModel_WithGenes_Reduced'],'s');
            end
            if strcmp(models{i},'Sobotta_Frontiers2017')==1  % just keep Setup_Core.m
                delete([target_path,filesep,models{i},filesep,'Setup_IL6_ExpDesign.m']);
                delete([target_path,filesep,models{i},filesep,'Setup_IL6_Full.m']);
            end
            if strcmp(models{i},'Swameye_PNAS2003')==1  % just keep Setup.m
                delete([target_path,filesep,models{i},filesep,'Setup2003.m']);                
            end
    end
end



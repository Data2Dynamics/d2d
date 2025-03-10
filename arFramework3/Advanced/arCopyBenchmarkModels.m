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
%               'lhsok' the models which can be fitted (20 benchmarks
%               except for Beer and Chen where LHS is pretty bad
% 
%               'fast10' 
%               the 10 models with least computation times
%
%               'new'
%               Three models that have been added to benchmark collection after
%               publication of Hass et al, 2018
%
%               'all' or '23'
%               The 20 original benchmark models plus the three new models
% 
%               The (starting characters of) model names can also be
%               specified, e.g.:
%               arCopyBenchmarkModels('Becker')
%               arCopyBenchmarkModels('Swa')
%
%               subset can also be a cell array of keywords and modelnames
%               e.g.: {'fast10', 'new', 'Bachmann'}.
% 
%   target_path     Default: current working directory (pwd)
% 
% Examples
% arCopyBenchmarkModels              % 20 default benchmark models
% arCopyBenchmarkModels('Bachmann')  % only the Bachmann model
% arCopyBenchmarkModels('Ra')        % only the Raia model (abbreviation Ra only fits for Raia)
% arCopyBenchmarkModels({'fast10', "Chen"})  % 10 fastest models plus Chen_MSB2009
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

function arCopyBenchmarkModels(subset, target_path)

if ~exist('target_path','var') || isempty(target_path)
    target_path = pwd;
end

if ~exist('subset','var') || isempty(subset)
    % set default
    subset = '20';

elseif iscell(subset)
    % recursive call
    cellfun(@(x) arCopyBenchmarkModels(x, target_path), subset);
    return

elseif isnumeric(subset)
    % convert to string
    subset = num2str(subset);
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


% models from benchmarking collection not yet included in Hass et al, 2019
new_models = {
              'Borghans_BiophysChem_1997', ...
              'Elowitz_Nature2000', ...
              'Sneyd_PNAS2002'};
all_models = sort([models20, new_models]);


switch lower(subset)
    case {'20',''} % default
        models = models20;
    case 'fast10' % the 10 fastest models
        models = {
            'Becker_Science2010',...
            'Boehm_JProteomeRes2014',...
            'Brannmark_JBC2010',...
            'Bruno_Carotines_JExpBio2016',...
            'Crauste_ImmuneCells_CellSystems2017',...
            'Fiedler_BMC2016',...
            'Fujita_SciSignal2010',...
            'Reelin_PONE2017',...
            'Schwen_InsulinMouseHepatocytes_PlosOne2014',...
            'Swameye_PNAS2003'};        
    case '5'
        models = {
            'Bachmann_MSB2011',...
            'Becker_Science2010',...
            'Beer_MolBiosyst2014',...
            'Raia_CancerResearch2011',...
            'Swameye_PNAS2003'};        
    case 'lhsok'
        models = {
            'Bachmann_MSB2011',...
            'Becker_Science2010',...
            'Boehm_JProteomeRes2014',...
            'Brannmark_JBC2010',...
            'Bruno_Carotines_JExpBio2016',...
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
    case 'new'
        models = new_models;
    case {'all', '23'}
        models = all_models;
    otherwise
        % find models that start with the string subset
        ind = startsWith(all_models, subset, IgnoreCase=true);
        models = all_models(ind);
        if isempty(models)
            error('Subset unknown');
        end
end

DoCopy(models, target_path);





function DoCopy(models, target_path)

source_path = fullfile(fileparts(which('arInit')), 'Examples');

for i=1:length(models)

    fprintf('Copying %s ...\n',models{i});
    copyfile(fullfile(source_path, models{i}), ...
             fullfile(target_path, models{i}));

    % model-specific operations which have to be done:
    switch models{i}
        case 'Becker_Science2010'
            % use only original model, not Setup_Regularization.m
            delete(fullfile(target_path, models{i}, 'Setup_Regularization.m'));

        case 'Boehm_JProteomeRes2014'
            % only keep Setup.m
            delete(fullfile(target_path, models{i}, 'Setup_FullModel_Boehm2014.m'));
        
        case 'Merkle_JAK2STAT5_PCB2016'
            % only keep SetupFinal.m
            delete(fullfile(target_path, models{i}, 'SetupCFUE.m'));
            delete(fullfile(target_path, models{i}, 'SetupComprehensive.m'));
            delete(fullfile(target_path, models{i}, 'SetupH838.m'));
            delete(fullfile(target_path, models{i}, 'SetupSens.m'));
        
        case 'Lucarelli_TGFb_2017'
            % just keep reduced model with genes
            rmdir(fullfile(target_path, models{i}, 'TGFb_ComplexModel_ModelSelection_L1'),'s');
            rmdir(fullfile(target_path, models{i}, 'TGFb_ComplexModel_withGenes'),'s');
            movefile(fullfile(target_path, models{i}, 'TGFb_ComplexModel_WithGenes_Reduced', '*'), ...
                     fullfile(target_path, models{i}));
            rmdir(fullfile(target_path, models{i}, 'TGFb_ComplexModel_WithGenes_Reduced'),'s');

        case 'Sobotta_Frontiers2017'
            % just keep Setup_Core.m
            delete(fullfile(target_path, models{i}, 'Setup_IL6_ExpDesign.m'));
            delete(fullfile(target_path, models{i}, 'Setup_IL6_Full.m'));

        case 'Swameye_PNAS2003'
            % just keep Setup.m
            delete(fullfile(target_path, models{i}, 'Setup2003.m'));
    end
    

end



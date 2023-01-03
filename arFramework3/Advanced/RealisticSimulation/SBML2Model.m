function SBML2Model     

%% Compile database model
if ~isempty(dir('*.xml'))
    Model = dir('*.xml');
    for i=1:length(Model)
        ModelName = Model(i).name(1:end-4);        
        if ~exist(['Models/' ModelName '.def'],'file')
            %arImportSBMLold(ModelName,'tend',100);
            arParseSBML(ModelName,'overwrite','tend',1000)
        end
    end
    arInit
    ar.config.checkForNegFluxes = false;
    for i=1:length(Model)
        arLoadModel(Model(i).name(1:end-4));
    end
    arCompileAll(true)
    arPlot
elseif exist('Setup','file')
    Setup
else
    error('Where is the Biomodel? Please help me to find it.')
end
arSave('Biomodel')
[~,ws]=fileparts(ar.config.savepath);
movefile(['Results' filesep ws],['Results' filesep 'Biomodel']);
fprintf('Biomodel workspace saved to file ./Results/Biomodel/workspace.mat \n');
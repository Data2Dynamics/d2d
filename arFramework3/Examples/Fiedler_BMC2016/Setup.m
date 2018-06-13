%% Compile model
arInit
arLoadModel('model_RafMekErk');
arLoadData('Ctrl',1);
arLoadData('Sorafenib_5_muM',1);
arLoadData('UO126_30_muM',1);
arCompileAll;

%% Set pre-equilibration
arSteadyState(1,arFindCondition(ar,'Ctrl'),1);

%% Parameter settings
arLoadPars('bestFit')


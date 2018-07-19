%% Compile model
arInit
arLoadModel('model_RafMekErk');
arLoadData('Ctrl',1);
arLoadData('Sorafenib_5_muM',1);
arLoadData('UO126_30_muM',1);
arCompileAll;

%% Parameter settings
arLoadPars('bestFit')

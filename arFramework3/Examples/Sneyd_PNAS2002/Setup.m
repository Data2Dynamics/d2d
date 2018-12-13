%% Compile model
arInit
arLoadModel('model_Sneyd_PNAS2002');
arLoadData('dataset_Sneyd_PNAS2002__fig4',1);
arLoadData('dataset_Sneyd_PNAS2002__fig5',1);
arCompileAll;

%% Load Parameters
arLoadPars('best_fit')

arPlot

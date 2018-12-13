%% loading model and a wild-type data-set where everything is observed
arInit;

%% Choose on of the three models:
arLoadModel('model1_dream6');
arLoadData('all_model1', 1, 'xls',false);      % everything observed

% arLoadModel('model2_dream6');
% arLoadData('all_model2', 1, 'xls',false);        % everything observed

% arLoadModel('model3_dream6');
% arLoadData('all_model3', 1, 'xls',false);      % everything observed

arCompileAll

ar.config.fiterrors = 0;

%% Loading of the original parameters and simulation of data:
loadGoldStandardParameter(ar.model(1).name) % no priors



% Load models & data
arInit

%Load full model or reduced
arLoadModel('model_expl4_full');
% arLoadModel('model_expl4_red');

%Load Data
arLoadData('data_expl4');
arCompileAll();


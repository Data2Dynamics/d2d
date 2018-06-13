% Load models & data
arInit

%Load full or reduced model
arLoadModel('model_expl3_full');
% arLoadModel('model_expl3_red');

%Load data set
arLoadData('data_expl3');
arCompileAll();
%Calibrate model, find global optimum

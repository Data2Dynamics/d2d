%load the d2d framework
arInit;

%load either the full or the reduced model
arLoadModel('model_expl1_full');
%arLoadModel('model_expl1_red');
%load the data set
arLoadData('data_expl1');

arCompileAll;
ar.config.useNewPlots=1;


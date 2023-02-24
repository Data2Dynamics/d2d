path_to_add = strrep(which('arInit.m'),'arInit.m',strrep('Examples\ToyModels\TransientFunction\TransientFunction_library','\',filesep));
addpath(path_to_add) % functions for fitting the transient model (e.g. fitting both signums)

%%
arInit
arLoadModel('model_template')
arLoadData('data_template')
arCompileAll

Initialize_FitTransient % setting bounds, identification of parameter-types from name

arFitLHS(10) % multi-start fitting
arFitTransient % single fit
arQplot('y')
arPlot


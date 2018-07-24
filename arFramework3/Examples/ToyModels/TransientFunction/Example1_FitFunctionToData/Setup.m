addpath('../TransientFunction_library/') % functions for fitting the transient model (e.g. fitting both signums)

%%
arInit
arLoadModel('TransientFunction')
arLoadData('ExampleData',1)
arCompileAll

Initialize_FitTransient % setting bounds, identification of parameter-types from name

arFitLHS(10) % multi-start fitting
arFitTransient % single fit
arQplot('y')
arPlot


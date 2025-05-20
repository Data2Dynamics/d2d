% Setup the model used for simulating only (no data included)

arInit;
arLoadModel('PHH_model_sim');
arLoadData('sim',1,'csv',true);
arLoadData('steadystate/steadystate_phh',1,'csv',true);
arCompileAll;
arFindInputs;
arSteadyState(1,arFindCondition(ar,'steady'),'all');

arDisableData(arFindData('steadystate','names'));
arDisableData(arFindData('sim','names'));

ar.model.qPlotYs=0*ar.model.qPlotYs;
ar.model.qPlotXs=[1,0];
save('Simulation_model_PHH');
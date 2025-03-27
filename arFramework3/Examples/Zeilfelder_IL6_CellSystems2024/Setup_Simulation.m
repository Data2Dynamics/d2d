arInit;
arLoadModel('combined_model');
arLoadData('sim',1,'csv',true);
arLoadData('steadystate/steadystate_dcf_w_smad',1,'csv',true);
arCompileAll;
arFindInputs;
arSteadyState(1,arFindCondition(ar,'steady'),'all');

arDisableData(arFindData('steadystate','names'));
arDisableData(arFindData('sim','names'));

ar.model.qPlotYs=0*ar.model.qPlotYs;
ar.model.qPlotXs=[1,0];
save('Simulation_model');
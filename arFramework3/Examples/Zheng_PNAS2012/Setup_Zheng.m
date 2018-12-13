% Load models & data
arInit;
arLoadModel('ModelHistoneZheng');
arLoadData('MethylationData_aggregated_NTKO');
arCompileAll;

arFindInputs;
arSteadyState(1,1,1,-1e7)
arLoadPars('bestFit')
arPrint
arPlot
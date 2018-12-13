arInit;
arLoadModel('detection_limit');
arLoadData('LoD');
arCompileAll(true);

ar.config.useNewPlots = 0;

ar.ub=ar.ub + 5;

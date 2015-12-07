arInit;
arLoadModel('Simplified_MAP_Kinase');
arLoadData('Simplified_MAP_Kinase');

arCompileAll;

ar.model(1).qPlotXs(:)=1;
ar.config.useNewPlots=1;
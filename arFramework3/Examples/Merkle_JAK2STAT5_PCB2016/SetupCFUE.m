% Load models & data

arInit;

arLoadModel('jak2_stat5_cfue');

arLoadData('CFUE/Long');
arLoadData('CFUE/Concentrations');
arLoadData('CFUE/RNA');
arLoadData('CFUE/ActD');
arLoadData('CFUE/Fine');
arLoadData('CFUE/SOCS3oe');
arLoadData('CFUE/SHP1oe');
arLoadData('CFUE/DoseResp_7min');
arLoadData('CFUE/DoseResp_30min');
arLoadData('CFUE/DoseResp_pSTAT5_10min_fine');
arLoadData('CFUE/DoseResp_CIS_90min');
arLoadData('CFUE/RNA_ActD_noEpo');

arCompileAll;

% parameters with prior information
[tmpmean,tmpsd] = arLogNtoN(4.1513, 4.1513*0.3); 
arSetPars('init_EpoRJAK2', tmpmean, [], [], [], [], 1, tmpmean, tmpsd);

% Fix RNA scale since not ID anyhow
arSetPars({'SOCS3RNAeqm','CISHRNAeqm'}, [0 0], [2 2]);

% load best fit parameter values
arLoadPars('bestFitCFUE');

arPlot;
arPrint;

ar.config.optim.Display = 'iter';
ar.model.qPlotYs(:) = 1;
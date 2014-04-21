% Load models & data

arInit;

arLoadModel('jak2_stat5_feedbacks');

arLoadData('CFUE_Long');
arLoadData('CFUE_Concentrations');
arLoadData('CFUE_RNA');
arLoadData('CFUE_ActD');
arLoadData('CFUE_Fine');
arLoadData('CFUE_CISoe');
arLoadData('CFUE_CISoe_pEpoR');
arLoadData('CFUE_SOCS3oe');
arLoadData('CFUE_SHP1oe');
arLoadData('CFUE_DoseResp_7min');
arLoadData('CFUE_DoseResp_30min');
arLoadData('CFUE_DoseResp_pSTAT5_10min_fine');
arLoadData('CFUE_DoseResp_CIS_90min');

arCompileAll;

% parameters with prior information
[tmpmean,tmpsd] = arLogNtoN(4.1513, 4.1513*0.3); 
arSetPars('init_EpoRJAK2', tmpmean, [], [], [], [], 1, tmpmean, tmpsd);

% Fix RNA scale since not ID anyhow
arSetPars({'SOCS3RNAEqc','CISRNAEqc'}, [0 0], [2 2]);

% load best fit parameter values
arLoadPars('BestFit');

arPlot;
arPrint;
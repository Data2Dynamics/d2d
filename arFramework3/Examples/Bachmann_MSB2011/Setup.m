% Load models & data

arInit;

arLoadModel('jak2_stat5_feedbacks', 1);

arLoadData('CFUE_Long', 1);
arLoadData('CFUE_Concentrations', 1);
arLoadData('CFUE_RNA', 1);
arLoadData('CFUE_ActD', 1);
arLoadData('CFUE_Fine', 1);
arLoadData('CFUE_CISoe', 1);
arLoadData('CFUE_CISoe_pEpoR', 1);
arLoadData('CFUE_SOCS3oe', 1);
arLoadData('CFUE_SHP1oe', 1);
arLoadData('CFUE_DoseResp_7min', 1);
arLoadData('CFUE_DoseResp_30min', 1);
arLoadData('CFUE_DoseResp_pSTAT5_10min_fine', 1);
arLoadData('CFUE_DoseResp_CIS_90min', 1);

arParseModel;
arWriteCFiles;
arLink;

% parameters with prior information
[tmpmean,tmpsd] = arLogNtoN(4.1513, 4.1513*0.3); 
arSetPars('init_EpoRJAK2', tmpmean, [], [], [], [], 1, tmpmean, tmpsd);

% Fix RNA scale since not ID anyhow
arSetPars({'SOCS3RNAEqc','CISRNAEqc'}, [0 0], [2 2]);

% load best fit parameter values
arLoadPars('BestFit');

arPlot;
arPrint;
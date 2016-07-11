% Load models & data

arInit;

arLoadModel('jak2_stat5_h838');

% H838
arLoadData('H838/IB_TC_pEpoR-pSTAT5');
arLoadData('H838/IB_TC_pJAK2-pEpoR');
arLoadData('H838/IB_TC_pJAK2-pSTAT5');
arLoadData('H838/IB_TC_pSTAT5');
arLoadData('H838/IB_TC_pJAK2-pEpoR-pSTAT5');
arLoadData('H838/MS_TC_pSTAT5');
arLoadData('H838/MS_DR_pSTAT5');
arLoadData('H838/qRTPCR_TC_CISHRNA');

% H838-EpoR
arLoadData('H838EpoR/IB_TC_pSTAT5-pEpoR');
arLoadData('H838EpoR/IB_TC_pJAK2-pEpoR');
arLoadData('H838EpoR/IB_TC_STAT5_1');
arLoadData('H838EpoR/IB_TC_STAT5_2');
arLoadData('H838EpoR/IB_DR_pJAK2-EpoR');
arLoadData('H838EpoR/IB_TC_SOCS3');
arLoadData('H838EpoR/IB_DR_SOCS3');
arLoadData('H838EpoR/MS_DR_pSTAT5');
arLoadData('H838EpoR/IB_DR_pJAK2-pEpoR-pSTAT5_1_epo');
arLoadData('H838EpoR/IB_DR_pJAK2-pEpoR-pSTAT5_2_epo');
arLoadData('H838EpoR/qRTPCR_TC_Cish-SOCS3');
arLoadData('H838EpoR/qRTPCR_ActD_noEpo');

% Both H838 and H838-EpoR on one blot
arLoadData('H838-H838EpoR_IB_DR_pJAK2-pEpoR');

arCompileAll;

% Fix RNA scale since not ID anyhow
arSetPars({'SOCS3RNAEqc','CISRNAEqc'}, [0 0], [2 2]);

% load best fit parameter values
arLoadPars('bestFitH838');

arPlot;
arPrint;

ar.config.optim.Display = 'iter';
ar.model.qPlotYs(:) = 1;
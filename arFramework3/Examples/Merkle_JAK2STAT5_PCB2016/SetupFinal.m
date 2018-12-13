% Load models & data

arInit;

arLoadModel('jak2_stat5_cfue');
arLoadModel('jak2_stat5_h838_l1_final');

% CFU-E
arLoadData('CFUE/Long',1);
arLoadData('CFUE/Concentrations',1);
arLoadData('CFUE/RNA',1);
arLoadData('CFUE/ActD',1);
arLoadData('CFUE/Fine',1);
arLoadData('CFUE/SOCS3oe',1);
arLoadData('CFUE/SHP1oe',1);
arLoadData('CFUE/DoseResp_7min',1);
arLoadData('CFUE/DoseResp_30min',1);
arLoadData('CFUE/DoseResp_pSTAT5_10min_fine',1);
arLoadData('CFUE/DoseResp_CIS_90min',1);
arLoadData('CFUE/RNA_ActD_noEpo',1);

% H838
arLoadData('H838_L1/IB_TC_pEpoR-pSTAT5',2);
arLoadData('H838_L1/IB_TC_pJAK2-pEpoR',2);
arLoadData('H838_L1/IB_TC_pJAK2-pSTAT5',2);
arLoadData('H838_L1/IB_TC_pSTAT5',2);
arLoadData('H838_L1/IB_TC_pJAK2-pEpoR-pSTAT5',2);
arLoadData('H838_L1/MS_TC_pSTAT5',2);
arLoadData('H838_L1/MS_DR_pSTAT5',2);
arLoadData('H838_L1/qRTPCR_TC_CISHRNA',2);

% H838-EpoR
arLoadData('H838EpoR_L1/IB_TC_pSTAT5-pEpoR',2);
arLoadData('H838EpoR_L1/IB_TC_pJAK2-pEpoR',2);
arLoadData('H838EpoR_L1/IB_TC_STAT5_1',2);
arLoadData('H838EpoR_L1/IB_TC_STAT5_2',2);
arLoadData('H838EpoR_L1/IB_DR_pJAK2-EpoR',2);
arLoadData('H838EpoR_L1/IB_TC_SOCS3',2);
arLoadData('H838EpoR_L1/IB_DR_SOCS3',2);
arLoadData('H838EpoR_L1/MS_DR_pSTAT5',2);
arLoadData('H838EpoR_L1/IB_DR_pJAK2-pEpoR-pSTAT5_1_epo',2);
arLoadData('H838EpoR_L1/IB_DR_pJAK2-pEpoR-pSTAT5_2_epo',2);
arLoadData('H838EpoR_L1/qRTPCR_TC_Cish-SOCS3',2);
arLoadData('H838EpoR_L1/qRTPCR_ActD_noEpo',2);

% Both H838 and H838-EpoR on one blot
arLoadData('H838-H838EpoR_L1_IB_DR_pJAK2-pEpoR',2);

% Validation
arLoadData('CFUE/RNA_validation',1);
arLoadData('H838EpoR_L1/qRTPCR_TC_Cish-SOCS3_long',2);

arCompileAll;
arFindInputs;
%ar.config.optim.Display = 'iter';
ar.model(1).qPlotYs(:) = 1;
ar.model(2).qPlotYs(:) = 1;

arLoadPars('reducedfinal');

% Not fit validation data
for im = 1:length(ar.model)
    for id = ar.model(im).plot(end).dLink
        ar.model(im).data(id).qFit(:) = 0;
    end
end

arPlot;
arPrint;
% Shaping the dose response of IL6 induced target gene transcription
global ar;

arInit;
modelname = 'il6_core';
arLoadModel(modelname);

svantje_app_calibration     = true;
svantje_app_validation      = true;
svantje_stattic             = true;
svantje_rux                 = true;
svantje_stattic_rux         = true;
bohl                        = true;
xiaoyun                     = true;

%% sebastian's data
if ( bohl )
    arLoadData('bohl/hep_2005_02_03_Cont40ng_T90min', 1, [], true);
    arLoadData('bohl/hep_2004_10_26_Cont_30min', 1, [], true);
    arLoadData('bohl/hep_2005_02_14_DoseResponse_T15min', 1, [], true);
    arLoadData('bohl/hep_2005_03_01_Pulse20min', 1, [], true);
    arLoadData('bohl/hep_2005_03_16_Cont40ng_23min', 1, [], true);
    arLoadData('bohl/hep_2005_07_26_Pulse5min', 1, [], true);
    arLoadData('bohl/hep_2005_11_22_DoseResponseUntil120ng', 1, [], true);
    arLoadData('bohl/hep_2009_10_14_qRTPCR', 1, [], true);
    arLoadData('bohl/mlc_per_cell_log', 1, [], true);
    arLoadData('bohl/hep_2007_07_10_ActD', 1, [], true);
    arLoadData('bohl/hep_2006_06_07_ActD', 1, [], true);
    arLoadData('bohl/hep_2006_02_09_FourPulses_ReceptorNucleus', 1, [], true);
    arLoadData('bohl/hep_2006_04_20_ThreePulses_ThreeDoses', 1, [], true);
    arLoadData('bohl/hep_2006_05_16_ThreeDoses', 1, [], true);
    arLoadData('bohl/hep_2006_11_15_FourPulses_120min', 1, [], true);
    arLoadData('bohl/hep_2007_04_04_DoseResponse_3TP', 1, [], true);
    arLoadData('bohl/hep_2007_07_17_DRTC', 1, [], true);
    arLoadData('bohl/hep_2009_06_23_DRTC', 1, [], true);
    arLoadData('bohl/hep_2009_07_07_DRTC', 1, [], true);
end

%% svantje's Stattic & Ruxolitinib data
if(svantje_stattic_rux)
    arLoadData('braun/hep_2013_11_04_Stattic_Ruxolitinib_Inhibitor', 1, [], true);
    arLoadData('braun/hep_2013_11_12_Stattic_Ruxolitinib_Inhibitor', 1, [], true);
    
    arLoadData('braun/app/hep_2014_04_28_IL6DR_Rux_1h_SOCS3', 1, [], true);
    
    % SOCS3
    arLoadData('braun/app/hep_2011_06_06_qPCR_140109_IL6DR_1h', 1, [], true);
    arLoadData('braun/app/hep_2013_10_14_qPCR_140109_IL6DR_1h', 1, [], true);
    arLoadData('braun/app/hep_2013_10_21_qPCR_140109_IL6DR_1h', 1, [], true);

    arLoadData('braun/app/hep_2012_02_14_qPCR_140224_IL6DR_Inh_1h_Socs3', 1, [], true);
    arLoadData('braun/app/hep_2012_04_10_qPCR_140224_IL6DR_Inh_1h_Socs3', 1, [], true);
    arLoadData('braun/app/hep_2013_10_14_qPCR_140224_IL6DR_Inh_1h_Socs3', 1, [], true);
    
    arLoadData('braun/app/hep_2014_04_22_qPCR_140526_IL6DR_Rux_6h_SOCS3', 1, [], true);
    arLoadData('braun/app/hep_2014_05_19_qPCR_140526_IL6DR_Rux_6h_SOCS3', 1, [], true);
    arLoadData('braun/app/hep_2014_04_22_qPCR_140604_IL6DR_Rux_24h_SOCS3', 1, [], true);
    arLoadData('braun/app/hep_2014_05_19_qPCR_140604_IL6DR_Rux_24h_SOCS3', 1, [], true);
end

%% svantje's data
arLoadData('braun/hep_2011_09_06_DRTC', 1, [], true);
arLoadData('braun/hep_2011_04_04_STAT3_pY_degree_replicates', 1, [], true);
arLoadData('braun/hep_2013_07_08_qPCR_Socs3_Stat3', 1, [], true);

%% svantje's Stattic data
if(svantje_stattic)
    arLoadData('braun/hep_2011_12_12_Stattic_Stat3_Inhibitor_DR', 1, [], true);
    arLoadData('braun/hep_2012_01_10_Stattic_Stat3_Inhibitor_DR', 1, [], true);
    arLoadData('braun/hep_2012_01_23_Stattic_Stat3_Inhibitor_TC', 1, [], true);
    d1 = length(ar.model(1).data)+1; 
	 arLoadData('braun/hep_2012_01_23_Stattic_Stat3_Inhibitor', 1, [], true);
	 d2 = length(ar.model(1).data); for d=d1:d2; arMakeMeanAndStd(1,d); end
    arLoadData('braun/hep_2012_02_14_Stattic_Stat3_Inhibitor_TC', 1, [], true);
    arLoadData('braun/hep_2012_04_10_Stattic_Stat3_Inhibitor_TC', 1, [], true);
    arLoadData('braun/hep_2012_05_02_Stattic_Stat3_Inhibitor_40', 1, [], true);
    arLoadData('braun/hep_2012_05_02_Stattic_Stat3_Inhibitor_60', 1, [], true);
end

%% svantje's Ruxolitinib data
if(svantje_rux)
    arLoadData('braun/hep_2013_06_03_Ruxolitinib', 1, [], true);
    arLoadData('braun/hep_2013_06_17_Ruxolitinib', 1, [], true);
    arLoadData('braun/hep_2013_09_30_Ruxolitinib_Inhibitor_TC', 1, [], true);
end

%% xiaoyun's data
if ( xiaoyun )
    arLoadData('xiaoyun/nucSTAT3_Ratio', 1, [], true);
end

%% equilibration condition
arLoadData('steadystate/steadystate',1,'csv',true)

%% compile
arCompileAll;
addpath('Helper');

% Tight equilibration tolerances
ar.config.eq_tol=1e-10;
ar.config.atol = 1e-8;
ar.config.rtol = 1e-8;
dataSets;

% Pre-equilibration settings (use a steady state condition to equilibrate the system)
arClearEvents;
arFindInputs;
arSteadyState(1,arFindCondition(ar,'steady'),'all');
arDisableData('steadystate_steadystate');

%% settings
ar.config.optimizer = 5;
ar.config.optim.TolX = 1e-8;
ar.config.optim.MaxIter = 1000;
ar.config.optim.Display = 'iter';
ar.config.maxsteps = 1000;
ar.config.fiterrors = 1;

arLoadPars('finalized_model');
%Set parameters to be estimated
ar.qFit(:) = 1;
arSetPars('scale_socs3_qpcr_braun_hep_2013_07_08_qPCR_Socs3_Stat3',[],0,1)
hill_loc = find(strcmp(ar.pLabel,'socs3rna_hill'));
arSetPars('socs3rna_hill',10.^ar.p(hill_loc),1,0,1,3)
arSetPars({'scale_pjak1_wb_bohl_hep_2009_07_07_DRTC','scale_pstat3_wb_bohl_hep_2009_07_07_DRTC','offset_pstat3_lumi_braun_hep_2012_05_02_Stattic_Stat3_Inhibitor_60','offset_pstat3_wb_bohl_hep_2005_02_14_DoseResponse_T15min','offset_pstat3_wb_bohl_hep_2005_07_26_Pulse5min','offset_pstat3_wb_braun_hep_2011_09_06_DRTC','offset_pstat3_wb_braun_hep_2012_05_02_Stattic_Stat3_Inhibitor_40','offset_pstat3_wb_braun_hep_2012_05_02_Stattic_Stat3_Inhibitor_60','offset_pstat3_wb_braun_hep_2013_06_17_Ruxolitinib'},zeros(1,9),zeros(1,9),zeros(1,9),ones(1,9)*-5,ones(1,9)*3);
arCalcMerit; arGetMerit;
%save(modelname, 'ar');

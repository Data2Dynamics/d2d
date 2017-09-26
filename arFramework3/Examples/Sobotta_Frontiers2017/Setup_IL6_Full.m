% Shaping the dose response of IL6 induced target gene transcription
global ar;

addpath('Helper');
arInit;
arLoadModel('il6_model');

svantje_stattic                 = true;
svantje_rux                     = true;
svantje_stattic_rux             = true;
svantje_app_calibration         = true;
svantje_app_validation          = true;
svantje_app_triple_validation   = true;
xiaoyun_validation              = true;
bohl                            = true;
xiaoyun                         = true;

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

%% svantje and xiaoyun's APP calibration data (use for fitting)
if(svantje_app_calibration)
    % CXCL10
    calibration{1} = 'braun/calibration/hep_2011_06_06_qPCR_140109_IL6DR_1h_CXCL10';
    calibration{2} = 'braun/calibration/hep_2013_10_14_qPCR_140109_IL6DR_1h_CXCL10';
    calibration{3} = 'braun/calibration/hep_2013_10_21_qPCR_140109_IL6DR_1h_CXCL10';
    
    calibration{4} = 'braun/calibration/hep_2014_01_13_IL6DR_6h_Replicate1';
    calibration{5} = 'braun/calibration/hep_2014_01_13_IL6DR_6h_Replicate2';
    calibration{6} = 'braun/calibration/hep_2014_01_13_IL6DR_6h_Replicate3';
    
    calibration{7} = 'braun/calibration/hep_2014_01_27_IL6DR_24h_Replicate1';
    calibration{8} = 'braun/calibration/hep_2014_01_27_IL6DR_24h_Replicate2';
    calibration{9} = 'braun/calibration/hep_2014_01_27_IL6DR_24h_Replicate3';
    
    calibration{10} = 'xiaoyun/calibration/APP_tc_replicate1';
    calibration{11} = 'xiaoyun/calibration/APP_tc_replicate2';
           
    arLoadData(calibration{1}, 1, [], true);
    arLoadData(calibration{2}, 1, [], true);
    arLoadData(calibration{3}, 1, [], true);
    arMergePlotMulti(1, { calibration{1}, calibration{2}, calibration{3}, }, ...
                        {'Replicate 1', 'Replicate 2', 'Replicate 3'}, 'APP_calibration_1h' );
    
    arLoadData(calibration{4}, 1, [], true);
    arLoadData(calibration{5}, 1, [], true);
    arLoadData(calibration{6}, 1, [], true);
    arMergePlotMulti(1, { calibration{4}, calibration{5}, calibration{6}, }, ...
                        {'Replicate 1', 'Replicate 2', 'Replicate 3'}, 'APP_calibration_6h' );
   
    arLoadData(calibration{7}, 1, [], true);
    arLoadData(calibration{8}, 1, [], true);
    arLoadData(calibration{9}, 1, [], true);
    arMergePlotMulti(1, { calibration{7}, calibration{8}, calibration{9}, }, ...
                        {'Replicate 1', 'Replicate 2', 'Replicate 3'}, 'APP_calibration_24h' );
                    
    arLoadData(calibration{10}, 1, 'csv', true);
    arLoadData(calibration{11}, 1, 'csv', true);
    
end

%% svantje's APP calibration data (use for validation; *not* fitting)
if(svantje_app_validation)
    validation{1} = 'braun/validation/Ruxolitinib_1h/hep_2013_10_14_qPCR_140224_IL6DR_Inh_1h_Cxcl10';
    validation{2} = 'braun/validation/Ruxolitinib_1h/hep_2014_04_28_IL6DR_Rux_1h_CXCL10';
    arLoadData( validation{1}, 1, [], true );
    arLoadData( validation{2}, 1, [], true );
    arMergePlot(1, validation{1}, validation{2}, 'Replicate 1', 'Replicate 2', 'APP_validation_1h' );
    
    validation{3} = 'braun/validation/Ruxolitinib_6h/hep_2014_04_22_qPCR_140526_IL6DR_Rux_6h_APP';
    validation{4} = 'braun/validation/Ruxolitinib_6h/hep_2014_05_19_qPCR_140526_IL6DR_Rux_6h_APP';
    arLoadData( validation{3}, 1, [], true );
    arLoadData( validation{4}, 1, [], true );
    arMergePlot(1, validation{3}, validation{4}, 'Replicate 1', 'Replicate 2', 'APP_validation_6h' );
    
    validation{5} = 'braun/validation/Ruxolitinib_24h/hep_2014_04_22_qPCR_140604_IL6DR_Rux_24h_APP';
    validation{6} = 'braun/validation/Ruxolitinib_24h/hep_2014_05_19_qPCR_140604_IL6DR_Rux_24h_APP';
    arLoadData( validation{5}, 1, [], true );
    arLoadData( validation{6}, 1, [], true );
    arMergePlot(1, validation{5}, validation{6}, 'Replicate 1', 'Replicate 2', 'APP_validation_24h' );

    validation{7} = 'braun/validation/Stattic_1h/hep_2012_02_14_qPCR_140224_IL6DR_Inh_1h_Cxcl10';
    validation{8} = 'braun/validation/Stattic_1h/hep_2012_04_10_qPCR_140224_IL6DR_Inh_1h_Cxcl10';
    arLoadData( validation{7}, 1, [], true );
    arLoadData( validation{8}, 1, [], true );
    arMergePlot(1, validation{7}, validation{8}, 'Replicate 1', 'Replicate 2', 'APP_validation_1h' );
end

%% svantje's APP triple validation data (use for validation; *not* fitting)
if(svantje_app_triple_validation)
    validation{9}  = 'braun/validation/Ruxolitinib_triple/2015-04-13_triple_treatment_replicate1';
    validation{10} = 'braun/validation/Ruxolitinib_triple/2015-04-13_triple_treatment_replicate2';
    validation{11} = 'braun/validation/Ruxolitinib_triple/2015-04-13_triple_treatment_replicate3';
    arLoadData( validation{9},  1, 'csv', true );
    arLoadData( validation{10}, 1, 'csv', true );
    arLoadData( validation{11}, 1, 'csv', true );
    arMergePlotMulti(1, {validation{9}, validation{10}, validation{11}}, {'Replicate 1', 'Replicate 2', 'Replicate 3'}, 'Triple_validation' );
end

if ( xiaoyun_validation )
    validation{12} = 'xiaoyun/validation/nucSTAT3_validation_20150428';
    validation{13} = 'xiaoyun/validation/nucSTAT3_validation_20150508';
    validation{14} = 'xiaoyun/validation/nucSTAT3_validation_20150528';
    arLoadData( validation{12}, 1, 'csv', true );
    arLoadData( validation{13}, 1, 'csv', true );
    arLoadData( validation{14}, 1, 'csv', true );
end

%% predictions
arLoadData('steadystate/steadystate',1,'csv',true)

%% compile
arCompileAll(true);

% Tight equilibration tolerances
ar.config.eq_tol=1e-10;
ar.config.atol = 1e-8;
ar.config.rtol = 1e-8;
dataSets;

% Pre-equilibration settings (use a steady state condition to equilibrate
% the system)
arClearEvents;
arFindInputs;
arSteadyState(1,arFindCondition(ar,'steady'),'all');

%% settings
ar.config.optimizer = 5;
ar.config.optim.TolX = 1e-8;
ar.config.optim.MaxIter = 1000;
ar.config.optim.Display = 'iter';
ar.config.maxsteps = 1000;
ar.config.fiterrors = 1;
arDisableData('steadystate_steadystate');

% Load estimated parameters
arLoadPars('finalized_model');

fn = 'final_IL6_model';
arCalcMerit; arGetMerit;
save(fn, 'ar');

arLoadPars('finalized_model');

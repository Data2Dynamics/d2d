clear all;
close all;
clc;

%% Compile model
arInit
arLoadModel('PKA_cycle');

arLoadData('JI09_150302_Drg345-343_CycNuc__4-ABnOH_and_ctrl',1);
arLoadData('JI09_150302_Drg345-343_CycNuc__Rp8_Br_cAMPS_pAB_and_ctrl',1);
arLoadData('JI09_150302_Drg345-343_CycNuc__Rp8_pCPT_cAMPS_pAB_and_ctrl',1);
arLoadData('JI09_150302_Drg345-343_CycNuc__Rp_cAMPS_pAB_and_ctrl',1);
arLoadData('JI09_150302_Drg345-343_CycNuc__4-ABnOH_and_Fsk',1);
arLoadData('JI09_150302_Drg345-343_CycNuc__Rp8_Br_cAMPS_pAB_and_Fsk',1);
arLoadData('JI09_150302_Drg345-343_CycNuc__Rp8_pCPT_cAMPS_pAB_and_Fsk',1);
arLoadData('JI09_150302_Drg345-343_CycNuc__Rp_cAMPS_pAB_and_Fsk',1);

arLoadData('JI09_150330_Drg350-348_CycNuc__ctrl',1);
arLoadData('JI09_150330_Drg350-348_CycNuc__IBMX100',1);
arLoadData('JI09_150330_Drg350-348_CycNuc__Fsk10',1);
arLoadData('JI09_150330_Drg350-348_CycNuc__IBMX100_and_Fsk10',1);
arLoadData('JI09_150330_Drg350-348_CycNuc__Sp8_Br_cAMPS_AM10',1);

arLoadData('JI09_150330_Drg353-351_CycNuc__4-ABnOH_and_ctrl',1);
arLoadData('JI09_150330_Drg353-351_CycNuc__Rp8_Br_cAMPS_pAB_and_ctrl',1);
arLoadData('JI09_150330_Drg353-351_CycNuc__Rp8_pCPT_cAMPS_pAB_and_ctrl',1);
arLoadData('JI09_150330_Drg353-351_CycNuc__Rp_cAMPS_pAB_and_ctrl',1);
arLoadData('JI09_150330_Drg353-351_CycNuc__4-ABnOH_and_Fsk',1);
arLoadData('JI09_150330_Drg353-351_CycNuc__Rp8_Br_cAMPS_pAB_and_Fsk',1);
arLoadData('JI09_150330_Drg353-351_CycNuc__Rp8_pCPT_cAMPS_pAB_and_Fsk',1);
arLoadData('JI09_150330_Drg353-351_CycNuc__Rp_cAMPS_pAB_and_Fsk',1);

arLoadData('JI09_151102_Drg421-418_Age__ctrl_and_ctrl',1);
arLoadData('JI09_151102_Drg421-418_Age__ctrl_and_Fsk1',1);
arLoadData('JI09_151102_Drg421-418_Age__ctrl_and_Fsk3',1);
arLoadData('JI09_151102_Drg421-418_Age__ctrl_and_Fsk10',1);
arLoadData('JI09_151102_Drg421-418_Age__H8925_and_ctrl',1);
arLoadData('JI09_151102_Drg421-418_Age__H894_and_Fsk3',1);
arLoadData('JI09_151102_Drg421-418_Age__H8910_and_Fsk3',1);
arLoadData('JI09_151102_Drg421-418_Age__H8925_and_Fsk3',1);

arLoadData('JI09_160201_Drg453-452_CycNuc__ctrl',1);
arLoadData('JI09_160201_Drg453-452_CycNuc__Fsk',1);
arLoadData('JI09_160201_Drg453-452_CycNuc__Sp8_Br_cAMPS_AM',1);

arLoadData('JI09_160126_Drg449-444_CycNuc__ctrl',1);
arLoadData('JI09_160126_Drg449-444_CycNuc__Fsk10_and_IBMX100',1);
arLoadData('JI09_160126_Drg449-444_CycNuc__Sp8_Br_cAMPS_AM10',1);

arCompileAll;

%% Numerical setting
ar.config.atol = 1e-6;
ar.config.rtol = 1e-6;
ar.config.maxsteps = 1e4;
ar.config.maxstepsize = 0.5;

ar.config.optim = optimset(...
    ar.config.optim,'TolFun',10^-8,'PrecondBandWidth',inf,...
    'display','iter-detailed');

%% Set pre-equilibration
arSteadyState(1,arFindCondition(ar,'JI09_150330_Drg350_348_CycNuc__ctrl'),1:length(ar.model.condition));
arSimu(true,true,true);

%% Parameter settings
% Global scaling factors and backgrounds
arSetPars('s_pRII_global',[],[],[],3,3.5);
arSetPars('s_pRII_Western',[],[],[],-1,5);
arSetPars('s_Calpha_global',[],[],[],3,3.5);
arSetPars('b_pRII_global'  ,[],[],[],-3,3,1,log10(100),0.1);
arSetPars('b_Calpha_global',[],[],[],-3,3,1,log10(100),0.1);

% Experiment specific scaling factors
arSetPars('s_pRII_JI09_150302_Drg345_343_CycNuc'  ,[],[],[],-3,3,1,0,0.1);
arSetPars('s_pRII_JI09_150330_Drg350_348_CycNuc'  ,[],[],[],-3,3,1,0,0.1);
arSetPars('s_pRII_JI09_150330_Drg353_351_CycNuc'  ,[],[],[],-3,3,1,0,0.1);
arSetPars('s_pRII_JI09_151102_Drg421_418_Age'     ,[],[],[],-3,3,1,0,0.1);
arSetPars('s_Calpha_JI09_160201_Drg453_452_CycNuc',0,0); % fixed as only one dataset

% Prior knowledge for different parameters
arSetPars('KD_IBMX'             ,[],[],[],[],[],1,log10(10)         ,0.2);
arSetPars('KD_Fsk'              ,[],[],[],[],[],1,log10(7)          ,0.2);
arSetPars('KD_cAMP'             ,[],[],[],[],[],1,log10(2.9)        ,0.2);
arSetPars('KD_H89'              ,[],[],[],[],[],1,log10(0.048)      ,0.2);
arSetPars('xi_pPDE'             ,[],[],[], 0, 3,1,log10(2.5)        ,3.0);
arSetPars('xi_pAC'              ,[],[],[],-3, 0,1,log10(1)          ,3.0);
arSetPars('xi_kf_RII_2__RII_C_2',[],[],[],-5, 0,1,log10(3.8e4/2.1e6),0.1);
arSetPars('xi_kf_RII_C_2__RII_2',[],[],[],-5, 5,1,log10(2.6e4/3.0e4),0.1);

% Weak regularization of differences in the faction of imported cAMP analogues
arSetPars('xi_i_Rp8_Br_cAMPS_pAB'  ,[],[],[],-5,0,1,0,3);
arSetPars('xi_i_Rp8_pCPT_cAMPS_pAB',[],[],[],-5,0,1,0,3);
arSetPars('xi_i_Rp_cAMPS_pAB'      ,[],[],[],-5,0,1,0,3);
arSetPars('xi_i_Sp8_Br_cAMPS_AM'   ,[],[],[],-5,0,1,0,3);

% Weak regularization of differences in binding rates of cAMP analogues
arSetPars('xi_b_Rp8_Br_cAMPS'      ,[],[],[],-5,5,1,0,3);
arSetPars('xi_b_Rp8_pCPT_cAMPS'    ,[],[],[],-5,5,1,0,3);
arSetPars('xi_b_Rp_cAMPS'          ,[],[],[],-5,5,1,0,3);
arSetPars('xi_b_Sp8_Br_cAMPS'      ,[],[],[],-5,5,1,0,3);

% Weak regularization of differences in KD values of cAMP analogues (in [mumol/l])
arSetPars('xi_KD_Rp8_Br_cAMPS'     ,[],[],[],-5,5,1,0,3);
arSetPars('xi_KD_Rp8_pCPT_cAMPS'   ,[],[],[],-5,5,1,0,3);
arSetPars('xi_KD_Rp_cAMPS'         ,[],[],[],-5,5,1,0,3);
arSetPars('xi_KD_Sp8_Br_cAMPS'     ,[],[],[],-5,5,1,0,3);

% Relative antibody binding of pRII amd C-aplha antibody
arSetPars('rel_open'   ,[],[],[],-5,0);
arSetPars('xi_rel_open',[],[],[],-5,0);

% Fixation of non-identifiable scales
arSetPars('PDE_total' ,0,0);
arSetPars('RII2_total',0,0);
arSetPars('AC_total'  ,0,0);

% Extension of bounds
arSetPars('kf_RII_2__RII_C_2',[],[],[],-5,5);

% Load best fit
arLoadPars('BestFit');

%% Model settings
% (here for model selected in Isensee et al., JCB, 2018)
model_options.partial_import = 0;
model_options.different_KD   = 1;
model_options.PDE_inhibition = 0;
model_options.AC_inhibition  = 0;

% Partial import
if ~model_options.partial_import
    arSetPars('xi_i_Rp8_Br_cAMPS_pAB',0,0);
    arSetPars('xi_i_Rp8_pCPT_cAMPS_pAB',0,0);
    arSetPars('xi_i_Rp_cAMPS_pAB',0,0);
    arSetPars('xi_i_Sp8_Br_cAMPS_AM',0,0);
end

% Different KD values
if ~model_options.different_KD
    arSetPars('xi_KD_Rp8_Br_cAMPS',0,0);
    arSetPars('xi_KD_Rp8_pCPT_cAMPS',0,0);
    arSetPars('xi_KD_Rp_cAMPS',0,0);
    arSetPars('xi_KD_Sp8_Br_cAMPS',0,0);
end

% Inhibition of PDE by the C-subunit
if ~model_options.PDE_inhibition
    arSetPars('xi_pPDE',0,0);
    arSetPars('kf_PDE_Csub',0,0,0);
    arSetPars('KD_PDE_Csub',0,0,0);
end

% Inhibition of AC by the C-subunit
if ~model_options.AC_inhibition
    arSetPars('xi_pAC',0,0);
    arSetPars('kp_AC',0,0,0);
    arSetPars('kdp_AC',0,0,0);
end

%% Parameter optimization
arFitLHS(1000);
arPlotChi2s;

%% Visualization
ar.model(1).qPlotYs = ones(size(ar.model(1).qPlotYs));

arMergePlotMulti(1,{'JI09_150302_Drg345-343_CycNuc__4-ABnOH_and_ctrl',...
                    'JI09_150302_Drg345-343_CycNuc__Rp8_Br_cAMPS_pAB_and_ctrl',...
                    'JI09_150302_Drg345-343_CycNuc__Rp8_pCPT_cAMPS_pAB_and_ctrl',...
                    'JI09_150302_Drg345-343_CycNuc__Rp_cAMPS_pAB_and_ctrl',...
                    'JI09_150302_Drg345-343_CycNuc__4-ABnOH_and_Fsk',...
                    'JI09_150302_Drg345-343_CycNuc__Rp8_Br_cAMPS_pAB_and_Fsk',...
                    'JI09_150302_Drg345-343_CycNuc__Rp8_pCPT_cAMPS_pAB_and_Fsk',...
                    'JI09_150302_Drg345-343_CycNuc__Rp_cAMPS_pAB_and_Fsk'},...
                   {'4-ABnOH/ctrl',...
                    'Rp8-Br-cAMPS-pAB(10)/ctrl',...
                    'Rp8-pCPT-cAMPS-pAB(10)/ctrl',...
                    'Rp-cAMPS-pAB(10)/ctrl',...
                    '4-ABnOH/Fsk(10)',...
                    'Rp8-Br-cAMPS-pAB(10)/Fsk(10)',...
                    'Rp8-pCPT-cAMPS-pAB(10)/Fsk(10)',...
                    'Rp-cAMPS-pAB(10)/Fsk(10)'},...
                   'JI09_150302_Drg345-343_CycNuc');

arMergePlotMulti(1,{'JI09_150330_Drg350_348_CycNuc__ctrl',...
                    'JI09_150330_Drg350_348_CycNuc__IBMX100',...
                    'JI09_150330_Drg350_348_CycNuc__Fsk10',...
                    'JI09_150330_Drg350_348_CycNuc__IBMX100_and_Fsk10',...
                    'JI09_150330_Drg350_348_CycNuc__Sp8_Br_cAMPS_AM10'},...
                   {'ctrl',...
                    'IBMX(100)',...
                    'Fsk(10)',...
                    'IBMX(100)/Fsk(10)',...
                    'Sp-8-Br-cAMPS-AM(10)'},...
                   'JI09_150330_Drg350_348_CycNuc');

arMergePlotMulti(1,{'JI09_150330_Drg353-351_CycNuc__4-ABnOH_and_ctrl',...
                    'JI09_150330_Drg353-351_CycNuc__Rp8_Br_cAMPS_pAB_and_ctrl',...
                    'JI09_150330_Drg353-351_CycNuc__Rp8_pCPT_cAMPS_pAB_and_ctrl',...
                    'JI09_150330_Drg353-351_CycNuc__Rp_cAMPS_pAB_and_ctrl',...
                    'JI09_150330_Drg353-351_CycNuc__4-ABnOH_and_Fsk',...
                    'JI09_150330_Drg353-351_CycNuc__Rp8_Br_cAMPS_pAB_and_Fsk',...
                    'JI09_150330_Drg353-351_CycNuc__Rp8_pCPT_cAMPS_pAB_and_Fsk',...
                    'JI09_150330_Drg353-351_CycNuc__Rp_cAMPS_pAB_and_Fsk'},...
                   {'4-ABnOH/ctrl',...
                    'Rp8-Br-cAMPS-pAB(10)/ctrl',...
                    'Rp8-pCPT-cAMPS-pAB(10)/ctrl',...
                    'Rp-cAMPS-pAB(10)/ctrl',...
                    '4-ABnOH/Fsk(10)',...
                    'Rp8-Br-cAMPS-pAB(10)/Fsk(10)',...
                    'Rp8-pCPT-cAMPS-pAB(10)/Fsk(10)',...
                    'Rp-cAMPS-pAB(10)/Fsk(10)'},...
                   'JI09_150330_Drg353-351_CycNuc');

arMergePlotMulti(1,{'JI09_151102_Drg421-418_Age__ctrl_and_ctrl',...
                    'JI09_151102_Drg421-418_Age__ctrl_and_Fsk1',...
                    'JI09_151102_Drg421-418_Age__ctrl_and_Fsk3',...
                    'JI09_151102_Drg421-418_Age__ctrl_and_Fsk10',...
                    'JI09_151102_Drg421-418_Age__H8925_and_ctrl',...
                    'JI09_151102_Drg421-418_Age__H894_and_Fsk3',...
                    'JI09_151102_Drg421-418_Age__H8910_and_Fsk3',...
                    'JI09_151102_Drg421-418_Age__H8925_and_Fsk3'},...
                   {'ctrl/ctrl',...
                    'ctrl/Fsk(1)',...
                    'ctrl/Fsk(3)',...
                    'ctrl/Fsk(10)',...
                    'H89(25)/ctrl',...
                    'H89(4)/Fsk(3)',...
                    'H89(10)/Fsk(3)',...
                    'H89(25)/Fsk(3)'},...
                   'JI09_151102_Drg421-418_Age');

arMergePlotMulti(1,{'JI09_160201_Drg453-452_CycNuc__ctrl',...
                    'JI09_160201_Drg453-452_CycNuc__Fsk',...
                    'JI09_160201_Drg453-452_CycNuc__Sp8_Br_cAMPS_AM'},...
                   {'ctrl',...
                    'Fsk(10)',...
                    'Sp8-Br-cAMPS-AM(10)'},...
                   'JI09_160201_Drg453-452_CycNuc');

arMergePlotMulti(1,{'JI09_160126_Drg449-444_CycNuc__ctrl',...
                    'JI09_160126_Drg449-444_CycNuc__Fsk10_and_IBMX100',...
                    'JI09_160126_Drg449-444_CycNuc__Sp8_Br_cAMPS_AM10'},...
                   {'ctrl',...
                    'Fsk(10)/IBMX(100)',...
                    'Sp8-Br-cAMPS-AM(10)'},...
                   'JI09_160126_Drg449-444_CycNuc');

arPlot

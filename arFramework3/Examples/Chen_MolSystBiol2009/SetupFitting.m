clear all;
close all;
clc;

%% Load Model
arInit
arLoadModel('erbb_signaling'); model_name = 'erbb_signaling';

%% Load Data
dsl = {};
dls_filename = {};

% experimentaldata1
dls_filename{end+1} = 'experimentaldata1';
dsl{end+1}.name = {'experimentaldata1'};

% experimentaldata2
dls_filename{end+1} = 'experimentaldata2';
dsl{end+1}.name = {'experimentaldata2'};

% experimentaldata3
dls_filename{end+1} = 'experimentaldata3';
dsl{end+1}.name = {'experimentaldata3'};

% experimentaldata4
dls_filename{end+1} = 'experimentaldata4';
dsl{end+1}.name = {'experimentaldata4'};

% Loop: Datasets and Conditions
warning off
for i = 1:length(dsl)
    for j = 1:length(dsl{i}.name)
        arLoadData(dsl{i}.name{j},1,'csv');
    end
end
warning on

% Compile model
warning off
arCompileAll;
warning on

% Show condition
arShowDataConditionStructure

%% Constraint parameters
ar.lb = ar.lb.*0-15;
ar.ub = ar.ub.*0+9;

% prior_scale__knowledge = 0.2;
% prior_scale__regularisation = 1;
% prior_scale__sparsity = 0.01;
% 
% % Total abundances
% arSetPars('ksyn_EGFR__MKN1_fm'  ,0,0);
% arSetPars('RAStotal__MKN1_fm'   ,0,0);
% arSetPars('MEKtotal__MKN1_fm'   ,0,0);
% arSetPars('MAPKtotal__MKN1_fm'  ,0,0);
% arSetPars('PI3Ktotal__HS746T_fm',0,0);
% arSetPars('MPI3Ktotal__MKN1_fm' ,0,0);
% arSetPars('AKTtotal__MKN1_fm'   ,0,0);
% 
% % Prior
% arSetPars('KD_EGFR_EGF'     ,[],[],[],-3,3,1,log10(2)             ,prior_scale__knowledge); % KleinMat2003, page 1
% arSetPars('KD_EGFR_CET'     ,[],[],[],-3,3,1,log10(0.39)          ,prior_scale__knowledge); % KimGor2008, page 2
% arSetPars('kdeg_EGFR'       ,[],[],[],-3,3,1,log10(log(2)/(12*60)),prior_scale__knowledge); % SorkinDue2010, page 19 (half-life 8-14 hr -> 12 hr used here)
% arSetPars('kimp_pEGFR_EGF_2',[],[],[],-3,3,1,log10(0.25)          ,prior_scale__knowledge); % SorkinDue2010, page 19 (k_e = 0.15 to 0.40 min)
% arSetPars('kdeg_MMET'       ,[],[],[],-3,3,1,log10(log(2)/(5*60)) ,prior_scale__knowledge); % http://www.ncbi.nlm.nih.gov/pubmed/2554238 (half-life 5 hr)
% arSetPars('ki_AKT'          ,[],[],[],-3,3,1,log10(1.5*60)        ,prior_scale__knowledge); % SchoeberlPac2009, kf83
% arSetPars('kdeg_EGFR'       ,[],[],[],-3,3,1,log10(24*60)         ,prior_scale__knowledge); % http://www.fasebj.org/content/early/2012/01/09/fj.11-197723.full.pdf
% arSetPars('kdeg_MMET'       ,[],[],[],-3,3,1,log10(24*60)         ,prior_scale__knowledge); % http://www.fasebj.org/content/early/2012/01/09/fj.11-197723.full.pdf
% 
% % Differences between receptor dimers
% arSetPars('xi_kdim_MMET'          ,[],[],[],-3,3,1,0,prior_scale__regularisation);
% arSetPars('xi_kimp_pMMET_2'       ,[],[],[],-3,3,1,0,prior_scale__regularisation);
% arSetPars('xi_kexp_pMMET_2_i'     ,[],[],[],-3,3,1,0,prior_scale__regularisation);
% arSetPars('xi_kdeg_pMMET_2_i'     ,[],[],[],-3,3,1,0,prior_scale__regularisation);
% arSetPars('xi_ka_RAS_pMMET_2'     ,[],[],[],-3,3,1,0,prior_scale__regularisation);
% arSetPars('xi_ka_PI3K_pMMET_2'    ,[],[],[],-3,3,1,0,prior_scale__regularisation);
% 
% arSetPars('xi_kdim_MMET_EGFR'     ,[],[],[],-3,3,1,0,prior_scale__regularisation);
% arSetPars('xi_kimp_pMMET_pEGFR'   ,[],[],[],-3,3,1,0,prior_scale__regularisation);
% arSetPars('xi_kexp_pMMET_pEGFR_i' ,[],[],[],-3,3,1,0,prior_scale__regularisation);
% arSetPars('xi_kdeg_pMMET_pEGFR_i' ,[],[],[],-3,3,1,0,prior_scale__regularisation);
% arSetPars('xi_ka_RAS_pMMET_pEGFR' ,[],[],[],-3,3,1,0,prior_scale__regularisation);
% arSetPars('xi_ka_PI3K_pMMET_pEGFR',[],[],[],-3,3,1,0,prior_scale__regularisation);
% 
% % Impact of mutation
% arSetPars('xi_ki_MPI3K',[],[],[],-3,0,1,-1,prior_scale__regularisation);
% 
% % Cell-line and medium differences
% prior_scale_2 = 0.2;
% % prior_type = 1; % normal distribution
% prior_type = 3; % Laplace distribution
% arSetPars('d_ksyn_MMET__fm_2_hm'      ,[],[],0,-3,3,prior_type,0,prior_scale__sparsity);
% arSetPars('d_ksyn_EGFR__fm_2_hm'      ,[],[],0,-3,3,prior_type,0,prior_scale__sparsity);
% arSetPars('d_RAStotal__fm_2_hm'       ,[],[],0,-3,3,prior_type,0,prior_scale__sparsity);
% arSetPars('d_MEKtotal__fm_2_hm'       ,[],[],0,-3,3,prior_type,0,prior_scale__sparsity);
% arSetPars('d_MAPKtotal__fm_2_hm'      ,[],[],0,-3,3,prior_type,0,prior_scale__sparsity);
% arSetPars('d_PI3Ktotal__fm_2_hm'      ,[],[],0,-3,3,prior_type,0,prior_scale__sparsity);
% arSetPars('d_MPI3Ktotal__fm_2_hm'     ,[],[],0,-3,3,prior_type,0,prior_scale__sparsity);
% arSetPars('d_AKTtotal__fm_2_hm'       ,[],[],0,-3,3,prior_type,0,prior_scale__sparsity);
% arSetPars('d_ksyn_EGFR__MKN1_2_HS746T',[],[],0,-3,3,prior_type,0,prior_scale__sparsity);
% arSetPars('d_RAStotal__MKN1_2_HS746T' ,[],[],0,-3,3,prior_type,0,prior_scale__sparsity);
% arSetPars('d_MEKtotal__MKN1_2_HS746T' ,[],[],0,-3,3,prior_type,0,prior_scale__sparsity);
% arSetPars('d_MAPKtotal__MKN1_2_HS746T',[],[],0,-3,3,prior_type,0,prior_scale__sparsity);
% arSetPars('d_AKTtotal__MKN1_2_HS746T' ,[],[],0,-3,3,prior_type,0,prior_scale__sparsity);
% 
% % Load parameters 
% arLoadPars('20170213T161410_model_14')
% 
% % Parameters for prediction
% % IDp1
% arSetPars('s_pEGFR_IDp1_HS746T_HM_1min__late_time_points_rep1' ,0,0);
% arSetPars('s_pEGFR_IDp1_HS746T_HM_3min__late_time_points_rep1' ,0,0);
% arSetPars('s_pEGFR_IDp1_HS746T_HM_15min__late_time_points_rep1',0,0);
% arSetPars('s_pEGFR_IDp1_HS746T_HM_4h__late_time_points_rep1'   ,0,0);
% arSetPars('s_pEGFR_IDp1_HS746T_HM_7h__late_time_points_rep1'   ,0,0);
% arSetPars('s_pEGFR_IDp1_HS746T_HM_22h__late_time_points_rep1'  ,0,0);
% arSetPars('s_pEGFR_IDp1_HS746T_HM_72h__late_time_points_rep1'  ,0,0);
% 
% % IDp2
% arSetPars('s_pEGFR_IDp2_MKN1_HM_1min__late_time_points_rep1' ,0,0);
% arSetPars('s_pEGFR_IDp2_MKN1_HM_3min__late_time_points_rep1' ,0,0);
% arSetPars('s_pEGFR_IDp2_MKN1_HM_15min__late_time_points_rep1',0,0);
% arSetPars('s_pEGFR_IDp2_MKN1_HM_4h__late_time_points_rep1'   ,0,0);
% arSetPars('s_pEGFR_IDp2_MKN1_HM_7h__late_time_points_rep1'   ,0,0);
% arSetPars('s_pEGFR_IDp2_MKN1_HM_22h__late_time_points_rep1'  ,0,0);
% arSetPars('s_pEGFR_IDp2_MKN1_HM_72h__late_time_points_rep1'  ,0,0);
% 
% % IDp3
% arSetPars('KD_METinh'      ,-10,0,1,-11,11);
% arSetPars('kdim_MMETinh'   , 10,0,1,-11,11);
% arSetPars('s_pEGFR_IDp3_HS746T_HM_1min__Met_inhibition_rep1'   ,0,0);
% arSetPars('s_pEGFR_IDp3_HS746T_HM_3min__Met_inhibition_rep1'   ,0,0);
% arSetPars('s_pEGFR_IDp3_HS746T_HM_15min__Met_inhibition_rep1'  ,0,0);
% arSetPars('s_pEGFR_IDp3_HS746T_HM_240min__Met_inhibition_rep1' ,0,0);
% 
% % IDp4
% arSetPars('KI50_siRNA_EGFR',0,0);
% arSetPars('s_pEGFR_IDp4_MKN1_HM_1min__EGFR_expression'   ,0,0);
% arSetPars('s_pEGFR_IDp4_MKN1_HM_3min__EGFR_expression'   ,0,0);
% arSetPars('s_pEGFR_IDp4_MKN1_HM_15min__EGFR_expression'  ,0,0);
% arSetPars('s_pEGFR_IDp4_MKN1_HM_240min__EGFR_expression' ,0,0);


%% Numerical settings
ar.config.atol = 1e-6;
ar.config.rtol = 1e-6;
ar.config.maxsteps = 1e4;

ar.config.optim = optimset(...
    ar.config.optim,'TolFun',10^-8,'PrecondBandWidth',inf,...
    'display','iter-detailed');

%% Set pre-equilibration
ic_MKN1_FM   = arFindCondition(ar,'MKN1_FM'  ,'MKN1'  ,1,'full_medium'  ,1,'EGF_level',0,'CET_level',0,'METinh_level',0,'siRNA_EGFR_level',0);
ic_MKN1_HM   = arFindCondition(ar,'MKN1_HM'  ,'MKN1'  ,1,'hunger_medium',1,'EGF_level',0,'CET_level',0,'METinh_level',0,'siRNA_EGFR_level',0);
ic_HS746T_FM = arFindCondition(ar,'HS746T_FM','HS746T',1,'full_medium'  ,1,'EGF_level',0,'CET_level',0,'METinh_level',0,'siRNA_EGFR_level',0);
ic_HS746T_HM = arFindCondition(ar,'HS746T_HM','HS746T',1,'hunger_medium',1,'EGF_level',0,'CET_level',0,'METinh_level',0,'siRNA_EGFR_level',0);
ic_MKN1_HM_EGFR_silencing_00_33 = arFindCondition(ar,'MKN1_HM'  ,'MKN1'  ,1,'hunger_medium',1,'EGF_level',0,'CET_level',0,'METinh_level',0,'siRNA_EGFR_level',0.33);
ic_MKN1_HM_EGFR_silencing_01_00 = arFindCondition(ar,'MKN1_HM'  ,'MKN1'  ,1,'hunger_medium',1,'EGF_level',0,'CET_level',0,'METinh_level',0,'siRNA_EGFR_level',1.0);
ic_MKN1_HM_EGFR_silencing_03_00 = arFindCondition(ar,'MKN1_HM'  ,'MKN1'  ,1,'hunger_medium',1,'EGF_level',0,'CET_level',0,'METinh_level',0,'siRNA_EGFR_level',3.0);
ic_MKN1_HM_EGFR_silencing_30_00 = arFindCondition(ar,'MKN1_HM'  ,'MKN1'  ,1,'hunger_medium',1,'EGF_level',0,'CET_level',0,'METinh_level',0,'siRNA_EGFR_level',30.0);

if ~isempty(ic_MKN1_FM)
    arSteadyState(1,ic_MKN1_FM  ,  arFindCondition(ar,'MKN1_FM'  ,'MKN1'  ,1,'full_medium'  ,1,'siRNA_EGFR_level',0));
end
if ~isempty(ic_MKN1_HM)
    arSteadyState(1,ic_MKN1_HM  ,  arFindCondition(ar,'MKN1_HM'  ,'MKN1'  ,1,'hunger_medium',1,'siRNA_EGFR_level',0));
end
if ~isempty(ic_HS746T_FM)
    arSteadyState(1,ic_HS746T_FM,  arFindCondition(ar,'HS746T_FM','HS746T',1,'full_medium'  ,1,'siRNA_EGFR_level',0));
end
if ~isempty(ic_HS746T_HM)
    arSteadyState(1,ic_HS746T_HM,  arFindCondition(ar,'HS746T_HM','HS746T',1,'hunger_medium',1,'siRNA_EGFR_level',0));
end
if ~isempty(ic_MKN1_HM_EGFR_silencing_00_33)
    arSteadyState(1,ic_MKN1_HM_EGFR_silencing_00_33,  arFindCondition(ar,'MKN1_HM'  ,'MKN1',1,'hunger_medium',1,'siRNA_EGFR_level',0.33));
end
if ~isempty(ic_MKN1_HM_EGFR_silencing_01_00)
    arSteadyState(1,ic_MKN1_HM_EGFR_silencing_01_00,  arFindCondition(ar,'MKN1_HM'  ,'MKN1',1,'hunger_medium',1,'siRNA_EGFR_level',1.0));
end
if ~isempty(ic_MKN1_HM_EGFR_silencing_03_00)
    arSteadyState(1,ic_MKN1_HM_EGFR_silencing_03_00,  arFindCondition(ar,'MKN1_HM'  ,'MKN1',1,'hunger_medium',1,'siRNA_EGFR_level',3.0));
end
if ~isempty(ic_MKN1_HM_EGFR_silencing_30_00)
    arSteadyState(1,ic_MKN1_HM_EGFR_silencing_30_00,  arFindCondition(ar,'MKN1_HM'  ,'MKN1',1,'hunger_medium',1,'siRNA_EGFR_level',30.0));
end

arSimu(true,true,true);


%% Parameters

% Fitting
arFit

% n_fit = 250;
n_fit = 100;
if n_fit >= 2
    arFitLHS(n_fit);
    arPlotChi2s
    
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 20 25])
    print('-depsc2','-r1200',['./Figures/' model_name '__optimizer_convergence']);
end
save all

% Print result
arPrint
%%
theta_true = [1.06386e-4 4.3309e-2 6.81902e4 3.54317 6.00588e1 1.81735e-4 4.98862e4 6.73816e-3 4.0749e-2 1.92391e-2 9.97194e-2 1.5543e-5 5.17473e-3 3.05684e-2 3.27962e-2 2.10189e-6 5.1794e-15 1.21498e-3 1.13102e-3 0.1 0.1 0.1];
theta_true = log10(theta_true);
arFits(theta_true)

%% Visualization
% ar.model(1).qPlotYs(end-5:end) = 1
arPlot
arPlotter


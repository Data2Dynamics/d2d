arInit;
name = 'control_SMAD';
arLoadModel(name);

% Omit Dorsomorphin data
cull = {'input_dorso', @(dorso)str2num(dorso)>0};

arLoadData( 'qpcr/qpcr_all',         1, 'csv', true, 'RemoveEmptyConds', 'RemoveConditions', cull ); % qPCR time courses
arLoadData( 'qpcr/qpcr_APAP',        1, 'csv', true, 'RemoveEmptyConds', 'RemoveConditions', cull ); % qPCR time courses in the presence of APAP
arLoadData( 'qpcr/qpcr_DR_BMP',      1, 'csv', true, 'RemoveConditions', cull ); % qPCR dose response w.r.t. BMP
arLoadData( 'wb_data/STAT3_TC',      1, 'csv', true, 'RemoveEmptyConds', 'RemoveConditions', cull ); % Western blot time courses
arLoadData( 'wb_data/SMAD_DR',       1, 'csv', true, 'RemoveConditions', cull ); % SMAD dose response data

%% Mock data file used tocontrol_SMAD_ID1_attempt2 equilibrate the system
arLoadData('steadystate/steadystate_dcf_w_smad',1,'csv',true);
arCompileAll(true);

%% Pre-equilibration settings (use the steady state condition to equilibrate the system)
arClearEvents;
arFindInputs;
arSteadyState(1,arFindCondition(ar,'steady'),'all');

arDisableData(arFindData('steadystate','names'));

%% Model requires tight equilibration tolerances
ar.config.eq_tol=1e-10;
ar.config.atol = 1e-10;
ar.config.rtol = 1e-10;

%% Optimization settings
ar.config.optimizer = 11;
ar.config.optim.TolX = 1e-8;
ar.config.optim.MaxIter = 1000;
ar.config.optim.Display = 'iter';
ar.config.maxsteps = 1000;
ar.config.maxsteps = 100*ar.config.maxsteps;

% Prevent standard deviations from being overestimated too much
ar.ub(arFindPar('sd_'))=0;

% qPCR scale is structurally non-identifiable => fix it
fixedScale = arFindPar('scale_smad6_qpcr_qpcr_qpcr_APAP_nExpID4', 'exact');
ar.qFit(fixedScale)=0;
ar.p(fixedScale)=0;

fixedScale = arFindPar('scale_ID3_qpcr_qpcr_qpcr_APAP_nExpID4', 'exact');
ar.qFit(fixedScale)=0;
ar.p(fixedScale)=0;

% Remove the extremely high dose of APAP since these cells were of questionable viability
removeAPAP50;

% Remove a few datapoints that were identified as outliers
removeBad;

% Set empirical bounds
ar.lb = -3*ones(size(ar.p));
ar.ub = 2*ones(size(ar.p));

ar.ub(arFindPar('scale_psmad')) = 4;
ar.lb(arFindPar('scale_psmad')) = -4;

ar.lb(arFindPar('scale_ID1_qpcr')) = -5;
ar.ub(arFindPar('scale_ID1_qpcr')) = 1;

ar.lb(arFindPar('smadinh_trans')) = -3;
ar.ub(arFindPar('smadinh_trans')) = 2;

ar.lb(arFindPar('DCF_mRNA_inh')) = -6;
ar.ub(arFindPar('DCF_mRNA_inh')) = 0;

ar.lb(arFindPar('Km_BMP')) = -1;
ar.ub(arFindPar('Km_BMP')) = 6;

ar.lb(arFindPar('sd_pSMAD')) = -1;
ar.ub(arFindPar('sd_pSMAD')) = 4;

ar.lb(arFindPar('smoothness_mRNA_loss')) = -1;
ar.ub(arFindPar('smoothness_mRNA_loss')) = 4;

% Remove Noggin data
removeNoggin;

arLink;

fixSelectedParameters;

% Better initial guess for optimization (the default doesn't evaluate properly).
ar.p = -1*ones(size(ar.p));

save(name, 'ar');
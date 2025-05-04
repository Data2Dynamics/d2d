% Set up the primary hepatocyte model.

arInit;
ar.config.checkForNegFluxes=2;

name = 'PHH_model';
arLoadModel(name);

%% Mock data file used to equilibrate the system
arLoadData('steadystate/steadystate_phh',1,'csv',true)

cull = {'input_dorso', @(dorso)str2num(dorso)>0, 'input_bmp', @(bmp)str2num(bmp)>0, 'input_bmp2', @(bmp2)str2num(bmp2)>0, 'is_bmp2', @(isbmp2)str2num(isbmp2)>0, 'time', @(t)(str2num(t)>240)};
arLoadData( 'PHH/2018_32_PHH_genes_raw',   1, 'xlsx', true, 'RemoveConditions', cull );
arLoadData( 'PHH/2018_32_PHH_pSTAT3',   1, 'xlsx', true, 'RemoveConditions', cull );
arLoadData( 'PHH/2018_42_PHH_MpC',   1, 'xlsx', true, 'RemoveConditions', cull );
arLoadData( 'PHH/2017_04_BSA',   1, 'xlsx', true, 'RemoveConditions', cull );
arLoadData( 'PHH/2017_04_PHH_Exp45', 1, 'xlsx', true, 'RemoveConditions', cull );
arLoadData( 'PHH/2018_49_SOCS3_IP_Exp137', 1, 'xlsx', true, 'RemoveConditions', cull );

arCompileAll;

%% Pre-equilibration settings (use the steady state condition to equilibrate the system)
arClearEvents;
arFindInputs;

tDefault = arFindCondition(ar,'steadystate_steadystate_phh', 'exact');
arSteadyState( 1, tDefault, 'all' );

arDisableData(arFindData('steadystate','names'));

%% Model requires tight equilibration tolerances
ar.config.eq_tol=1e-10;
ar.config.atol = 1e-8;
ar.config.rtol = 1e-8;

%% Optimization settings
ar.config.optimizer = 11;
ar.config.optim.TolX = 1e-10;
ar.config.optim.MaxIter = 1000;
ar.config.optim.Display = 'iter';
ar.config.maxsteps = 1000;
ar.config.maxsteps = 100*ar.config.maxsteps;

ar.p(arFindPar('input_phh'))=0;
ar.qFit(arFindPar('input_phh'))=0;
ar.qLog10(arFindPar('input_phh'))=0;

save(name, 'ar');
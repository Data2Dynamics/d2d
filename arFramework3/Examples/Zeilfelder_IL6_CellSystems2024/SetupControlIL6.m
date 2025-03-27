arInit;
ar.config.checkForNegFluxes=2;
arLoadModel('control_IL6');

% We omit the following data
% 1. Data acquired in the prescence of diclofenac.
% 2. Data acquired in the prescence of acetaminophen.
% 3. Data acquired in the prescence of BMP2.
cull = {'input_apap', @(apap)str2num(apap)>0, 'input_dcf', @(dcf)str2num(dcf)>0, 'input_bmp', @(bmp)str2num(bmp)>0, 'input_bmp2', @(bmp2)str2num(bmp2)>0, 'is_bmp2', @(isbmp2)str2num(isbmp2)>0, 'time', @(t)(str2num(t)>240)};

arLoadData( 'wb_data/STAT3_DR_IL6',             1, 'csv', true, 'RemoveConditions', cull );         % Western blot dose response w.r.t. IL6
arLoadData( 'wb_data/STAT3_TC',                 1, 'csv', true, 'RemoveConditions', cull );         % Western blot time courses
arLoadData( 'wb_data/STAT3_TC_APAP',            1, 'csv', true, 'RemoveConditions', cull );         % Western blot time courses in the presence of APAP
arLoadData( 'qpcr/qpcr_APAP',                   1, 'csv', true, 'RemoveConditions', cull );         % qPCR time courses in the presence of APAP
arLoadData( 'qpcr/qpcr_DR_IL6',                 1, 'csv', true, 'RemoveConditions', cull );         % qPCR dose responses w.r.t. IL6
arLoadData( 'qpcr/qpcr_all',                    1, 'csv', true, 'RemoveConditions', cull );         % qPCR time courses
arLoadData( 'SOCS3_wb/Anja_2018_13_SOCS3',      1, 'csv', true, 'RemoveConditions', cull );         % SOCS3 protein data
arLoadData( 'SOCS3_wb/Anja_2019_03_SOCS3',      1, 'csv', true, 'RemoveConditions', cull );         % SOCS3 protein data
arLoadData( 'BSA/protein_conc',                 1, 'csv', true );                                   % Protein amounts
arLoadData( 'gp130/gp130_degradation',          1, 'csv', true );                                   % gp130 degradation experiment
arLoadData( 'mlc/mlc',                          1, 'xls', true, 'RemoveConditions', cull );         % Molecules per cell measurement

%% Mock data file used to equilibrate the system
arLoadData('steadystate/steadystate_dcf',1,'csv',true)
arCompileAll(true);

%% Pre-equilibration settings (use the steady state condition to equilibrate the system)
arClearEvents;
arFindInputs;
arSteadyState(1,arFindCondition(ar,'steady'),'all');

% Don't include the steady state condition in the fit.
arDisableData(arFindData('steadystate','names'));

%% Model requires tight equilibration tolerances
ar.config.eq_tol=1e-10;
ar.config.atol = 1e-8;
ar.config.rtol = 1e-8;

%% Optimization settings
ar.config.optimizer = 11;
ar.config.optim.TolX = 1e-8;
ar.config.optim.MaxIter = 1000;
ar.config.optim.Display = 'iter';
ar.config.maxsteps = 1000;
ar.config.maxsteps = 100*ar.config.maxsteps;

% The hill coefficient must be more strictly bounded to avoid problems when
% optimizing.
ar.lb(arFindPar('n_STAT3'))=-1;
ar.ub(arFindPar('n_STAT3'))=log10(4);

% Fix specific parameters.
fixSelectedParameters;

%% Save model after Setup
save('control_IL6', 'ar');
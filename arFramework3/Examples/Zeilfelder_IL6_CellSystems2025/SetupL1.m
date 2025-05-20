arInit;
ar.config.checkForNegFluxes=2;

name = 'L1_model';

arLoadModel(name);

cull = {'input_dorso', @(dorso)str2num(dorso)>0, 'input_bmp', @(bmp)str2num(bmp)>0, 'input_bmp2', @(bmp2)str2num(bmp2)>0, 'is_bmp2', @(isbmp2)str2num(isbmp2)>0, 'time', @(t)(str2num(t)>240)};
     
arLoadData( 'wb_data/STAT3_DR_IL6',                 1, 'csv', true, 'RemoveConditions', cull );  %X
arLoadData( 'wb_data/STAT3_TC',                     1, 'csv', true, 'RemoveConditions', cull );  %X
arLoadData( 'wb_data/STAT3_TC_APAP',                1, 'csv', true, 'RemoveConditions', cull );  %X

arLoadData( 'qpcr/qpcr_APAP',                       1, 'csv', true, 'RemoveConditions', cull );
arLoadData( 'qpcr/qpcr_DR_IL6',                     1, 'csv', true, 'RemoveConditions', cull ); %X
arLoadData( 'qpcr/qpcr_all',                        1, 'csv', true, 'RemoveConditions', cull ); %X

arLoadData( 'SOCS3_wb/Anja_2018_13_SOCS3',          1,  'csv', true, 'RemoveConditions', cull );  %X
arLoadData( 'SOCS3_wb/Anja_2019_03_SOCS3',          1,  'csv', true, 'RemoveConditions', cull );
arLoadData( 'BSA/protein_conc',                     1,  'csv', true ); %X
arLoadData( 'gp130/gp130_degradation',              1,  'csv', true ); %X
arLoadData( 'mlc/mlc',                              1,  'xls',true);

%% Mock data file used to equilibrate the system
arLoadData('steadystate/steadystate_dcf',1,'csv',true)

arCompileAll;

%% Pre-equilibration settings (use the steady state condition to equilibrate the system)
arClearEvents;
arFindInputs;
arSteadyState(1,arFindCondition(ar,'steady'),'all');

% Don't include the steady state data in the fit
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

% Remove the extremely high dose of APAP since these cells were of questionable viability
% removeAPAP50;

% Remove a few datapoints that were identified as outliers
removeBad;

% Remove unused parameters
fixSelectedParameters;

% Relink the model
arLink;

save(name, 'ar');